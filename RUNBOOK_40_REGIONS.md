# Runbook — passage aux 40 régions, benchmark, et mise à jour des chiffres

Procédure pas à pas. Chaque étape est reprenable : rien n'est recalculé si le
résultat est déjà sur le disque, sauf `--force`.

Tous les scripts vivent dans `test/` et se lancent depuis ce répertoire.

---

## Prérequis

1. **Les deux dépôts côte à côte** (ou `HEALPIX_RESAMPLE_ROOT` pointant vers le
   second) :

   ```
   Projects/
     healpix-compress/
     healpix-resample/
   ```

   ```bash
   # si le checkout resample est ailleurs
   export HEALPIX_RESAMPLE_ROOT=/chemin/vers/healpix-resample
   ```

2. **Environnement** — la pile d'évaluation est désormais déclarée dans
   `pyproject.toml` sous l'extra `paper` :

   ```bash
   pip install -e ".[paper]"
   # + les codecs binaires si absents : zfpy, pysz, imagecodecs
   ```

   `healpix-resample` doit être importable : c'est lui qui fournit
   `CloughTocherResampler` (peuplement des tuiles) et `KNeighborsResampler.invert`
   (référence de reconstruction commune).

3. **Réseau** — seule l'étape 2 en a besoin (lecture des COG earth-search).
   Elle lit un bucket public en anonyme ; si `AWS_ACCESS_KEY_ID` traîne dans
   l'environnement, GDAL signe la requête et échoue : `unset` avant de lancer.

---

## Ce que font les nouveaux fichiers

| Fichier | Rôle |
|---|---|
| `paper_config.py` | Toutes les constantes géométriques et de protocole, en un seul endroit |
| `region_sites.py` | Lit les 40 régions dans le `data_manifest.csv` de healpix-resample |
| `build_region_windows.py` | Ré-extrait des fenêtres 1512² aux mêmes ancres, depuis les mêmes produits S2 |
| `build_tiles.py` | Tuile HEALPix + carré UTM équivalent-aire, avec assertions dures |
| `codec_lib.py` | Adaptateurs codecs et balayages — source unique (fin des 3 copies) |
| `run_codec_bench.py` | Balayage débit-distorsion, reprenable, deux grilles |
| `common_reference.py` | Opérateur direct HEALPix → pixels UTM (Limitation 5) |
| `paper_stats.py` | Interpolation log-log, codec oracle, Wilcoxon, Holm, Spearman |
| `update_paper_numbers.py` | Recalcule tous les chiffres et réécrit les tables LaTeX |
| `window_geometry.py` | Calcule, sans aucune image, la fenêtre que chaque ancre exige et où la centrer |
| `check_install.py` | Vérifie que la chaîne est complète et cohérente sur cette machine |
| `make_figures.py` | Régénère les figures du papier |

Les notebooks existants restent utilisables pour l'exploration ; ils ne sont
plus la source de vérité des chiffres.

---

## Étape 0 — vérifications, sans rien télécharger

```bash
cd test
python check_install.py         # tous les modules présents et importables ?
python paper_config.py          # 637 px, boîte englobante 901 px, fenêtre sûre 1801 px
python region_sites.py          # les 40 régions, classes, latitudes, produits
```

La chaîne est répartie sur une douzaine de modules voisins : un fichier manquant
après un transfert entre machines ne se voit qu'au milieu d'un run, sous la forme
d'un `ModuleNotFoundError` opaque. `check_install.py` les vérifie tous d'un coup,
avec taille et empreinte, et `python check_install.py --sha` donne des lignes
directement comparables entre deux machines :

```bash
# sur la machine de référence, puis sur l'autre, et on diffe
python check_install.py --sha > /tmp/here.txt
```

Le transfert fiable reste le dépôt git : ces fichiers sont neufs et donc non
suivis tant qu'ils n'ont pas été committés.

```bash
git add test/*.py *.md && git commit -m "paper pipeline" && git push
```

Attendu : `40 regions | classes {agriculture: 10, forest: 10, urban: 10, water: 10} |
HEALPix regime {equatorial: 22, cap: 18} | lat -43.5..64.1 | pinned products 36/40`.

Les 4 régions sans produit épinglé — black_forest, finland_boreal, lake_constance,
lake_geneva — sont un cas connu, et le script le signale de lui-même :
`healpix-resample` a échoué exactement sur les mêmes quatre, après 9 positions de
lattice chacune (voir `notebooks/tables/real_groundtruth_multiregion_products.csv`,
4 lignes `failed`), et s'est rabattu sur le catalogue EOPF. Voir la section
« Les quatre régions sans COG public » plus bas.

**Régression sur l'existant** — à faire avant tout le reste, c'est ce qui prouve
que la nouvelle chaîne reproduit le papier actuel :

```bash
python update_paper_numbers.py --no-write-tex
cat results/paper_numbers_self.md
```

Sur les 7 scènes actuelles, ce script redonne exactement les chiffres publiés :
ratio brut 1.988, ratios à RMSE cible 1.588 / 1.467 / 1.239 / 0.801 / 0.910,
ratios à budget 13.445 / 4.167 / 1.315 / 0.856, γ★ 2.06 / 1.91 / 1.61 / 1.04 / 1.18.
Deux écarts assumés, tous deux documentés dans le script :

* les p-values sont maintenant calculées sur `log(ratio)` et non sur `ratio − 1`
  (un ratio et son inverse pèsent alors pareil) ; l'ancienne convention reste
  stockée sous `p_linear` dans le JSON ;
* **la Table VI est fausse dans le papier actuel.** L'audit (`assert` sur le
  décompte) donne 65 points de fonctionnement et : jpeg2000\_lossy 26 victoires
  (rang moyen 1.56), jpegxl\_lossy 22 (mais seulement 28 points évaluables),
  sz3 8 (rang moyen 2.57), zfp 6, dwt\_zstd 3. Le papier annonce « SZ3 gagne
  45/51, rang moyen 1.47 » : vérification faite point par point sur les CSV,
  JPEG 2000 est plus petit que SZ3 à quasiment chaque cible (par ex. urban à
  RMSE 1e-4 : 925 ko contre 948 ko). La phrase « SZ3 is the strongest
  general-purpose codec across most of the range » est à retirer, pas seulement
  à corriger sur le dénominateur.

---

## Étape 1 — acquisition des fenêtres 1512² (réseau)

Toujours commencer par deux régions :

```bash
python build_region_windows.py --regions paris amazon
python build_region_windows.py --check
```

Puis les 40 :

```bash
python build_region_windows.py            # ~10-30 s par région
python build_region_windows.py --check
```

Sorties : `test/data/regions/<scene>_data.zarr` + `acquisition_log.csv`
(produit retenu, épinglé ou résolu, date, nuages, fraction de pixels valides).

Points à surveiller dans le journal :

* `n_rejected_candidates > 0` : des produits ont été écartés avant d'en trouver un
  lisible. La colonne `rejected_candidates` dit pourquoi.
* `collection` : si elle vaut autre chose que `sentinel-2-l2a`, la région vient
  d'une Collection-1 — à mentionner dans le papier (ligne de traitement
  différente de celle des 36 autres).

* `finite_fraction < 1` : la fenêtre touche un bord de granule ou du no-data.
  L'étape 2 rejettera la région ; c'est voulu.
* `product_resolution = pinned-missing->resolved` : le produit du manifest
  resample n'est plus servi par earth-search. À mentionner dans le papier si
  cela concerne plus d'une poignée de régions.

La fenêtre est centrée sur la **même ancre** que le patch 256² publié : ce dernier
en est exactement le carré central, ce qui rend les deux papiers directement
superposables.

---

## Étape 2 — tuiles HEALPix + crops UTM

```bash
python build_tiles.py --regions paris amazon      # d'abord deux
python build_tiles.py                             # les 40
```

Sorties : `test/data/tiles/<scene>.npz` + `tiles_index.csv`.

Le script **rejette** une région (au lieu de la rafistoler) si :

* aucune tuile parente niveau 10 n'est entièrement dans la fenêtre avec halo ;
* la fraction de no-data **à l'intérieur de l'empreinte de la tuile** dépasse
  `--max-nodata` (défaut 0,1 %) ;
* le carré UTM 637² ne tient pas dans la fenêtre, ou contient du no-data ;
* le resampler ne rend pas exactement les 1 048 576 cellules demandées.

C'est le garde-fou direct contre l'effet que vous aviez signalé : une tuile
partiellement remplie se compresserait presque gratuitement et gonflerait le côté
HEALPix.

Deux réglages de sélection, tous deux déterministes et indépendants de la mesure :

* **Choix de la tuile parente.** Les tuiles candidates sont classées par distance
  entre leur centre et le centre de la fenêtre, pas par nombre de pixels couverts,
  et **toutes** sont évaluées. Une fenêtre de 1512 px couvre ~2,4 tuiles par axe :
  la tuile qui contient le plus de pixels est souvent une tuile qui dépasse d'un
  côté, alors qu'une voisine tient entièrement.
* **Traitement du no-data.** Les échantillons no-data sont retirés de la
  triangulation plutôt que passés en NaN (un NaN empoisonne tous les triangles qui
  le touchent). Sous le plafond `--max-nodata`, une poignée de pixels isolés est
  donc interpolée ; au-dessus, la région est rejetée. `tiles_index.csv` enregistre
  `nodata_in_tile_fraction` et `nodata_samples_dropped` pour chaque région : à
  reporter dans le papier (nombre de régions concernées et fraction maximale).

Coût : la triangulation de Delaunay porte sur la sous-fenêtre (~950² points, pas
1512²), soit quelques minutes par région. `--halo` (défaut 32) élargit la marge si
une région passe juste.

Contrôle d'échelle optionnel, sur les 7 scènes historiques :

```bash
python build_tiles.py --legacy --fill-nan-mean
```

---

## Étape 3 — balayage débit-distorsion

Test à blanc, puis la totalité :

```bash
python run_codec_bench.py --limit 2 --grid both
python run_codec_bench.py --grid both
```

Sorties : `results/bench_parts/<grid>__<scene>.csv`, puis assemblage automatique
en `results/hpx_codec_bench.csv` et `results/utm_codec_bench.csv` (mêmes colonnes
qu'avant, plus `grid`, `scene_class`, `lat`, `abs_lat`, `hpx_regime`,
`raw_bytes`, `rmse_common`).

Interrompu ? Relancer la même commande : les régions déjà écrites sont sautées.
Pour ne réassembler que les CSV : `python run_codec_bench.py --assemble`.

**Lire le compte de points.** Chaque grille affiche désormais `N/M sweep points`
et écrit `results/bench_parts/<grille>__<région>.skips.json` détaillant chaque
configuration écartée et pourquoi. Un `WARNING` s'affiche si un codec ne produit
aucun point sur une grille : c'est ce qui transforme silencieusement la
comparaison « codec oracle » en « quel codec était évaluable » plutôt que
« lequel était meilleur ». À vérifier avant de faire tourner les 40.

**Formes décodées.** Un codec à transformée peut reconstruire une image plus
grande que l'entrée : PyWavelets, sur le crop UTM de 637×637 (impair), rend
638×638, parce que les bandes de détail de niveau 1 font ceil(637/2) = 319. Cela
n'apparaît jamais sur la tuile HEALPix, dont le côté est une puissance de deux.
La reconstruction est maintenant recadrée à la taille d'entrée, et
`decoded_shape` / `cropped_to` sont écrits dans le `config_json` : les octets
encodés incluent toujours le remplissage, seule la comparaison est recadrée.

**Versions des bibliothèques.** `results/bench_environment.json` enregistre les
versions de numpy, PyWavelets, zfpy, pysz, imagecodecs, Pillow, zstandard, etc.
Ce n'est pas cosmétique : le comportement d'un codec dépend de sa version, et le
CSV publié en juillet contient 49 points `dwt_zstd` sur des crops 637×637 que la
version actuelle de PyWavelets ne peut pas produire — l'environnement a donc
changé sous les chiffres du papier. Cela ferme l'item « codec library and version
not reported » des Limitations.

Le balayage JPEG XL a été étendu : l'échelle `level` sature vers RMSE 1e-3 sur du
B04, donc JPEG XL n'était **pas évaluable** aux cibles serrées et n'y gagnait ni
n'y perdait par mérite. Des points `distance` (2.0 → 0.01) ont été ajoutés pour
couvrir le haut de gamme. Vérifier après le premier run que `jpegxl_lossy` a bien
autant de points évaluables que les autres — c'est exactement la réserve n° 4 de
`review_soumission_tgrs.md`.

---

## Étape 4 — référence de reconstruction commune

C'est le poste qui peut déplacer les chiffres, donc à faire avant toute
réécriture des Résultats.

```bash
python run_codec_bench.py --grid hpx --common-reference --force
```

Ajoute la colonne `rmse_common` (HEALPix décodé → projeté sur les pixels UTM via
`KNeighborsResampler.invert`, comparé à la même observation que le côté UTM) et
écrit `results/common_reference_floor.csv` : le **plancher d'erreur de
représentation**, c'est-à-dire la RMSE de la tuile *non compressée* une fois
reprojetée. Aucune configuration de codec ne peut passer sous ce plancher ; toute
cible RMSE inférieure devient inatteignable côté HEALPix, et le papier doit le
dire au lieu de rapporter un `n` réduit sans explication.

Asymétrie de modèle à assumer par écrit : les tuiles sont peuplées en
Clough–Tocher (opérateur de référence retenu) alors que la projection retour
utilise l'opérateur d'observation gaussien. La tuile n'est donc pas l'inverse de
l'opérateur qui prédit, et le plancher est plus haut qu'il ne le serait avec
`PSFResampler` des deux côtés. Contrôle de sensibilité, si vous voulez le
chiffrer :

```bash
python run_codec_bench.py --grid hpx --common-reference --cr-mode psf --force
```

---

## Étape 5 — chiffres et figures

```bash
python update_paper_numbers.py                    # métrique interne (comme aujourd'hui)
python update_paper_numbers.py --metric common    # métrique commune
python make_figures.py
```

Écrit dans `../tex/tables/` : `table_storage_accounting`, `table_storage_at_constant_rmse`,
`table_error_at_constant_size`, `table_gamma_star`, `table_codec_ranking`,
`table_benchmark_regions`, `table_latitude_effect`, `table_fixed_codec`.

Écrit dans `../tex/figures/` : `fig_ratio_boxplots`, `fig_utm_crop_footprints`,
`fig_region_map`, `fig_latitude_effect`.

Et surtout `results/paper_numbers_{self,common}.md` : le digest à garder ouvert
pendant la rédaction (tous les ratios, p bruts, p de Holm, γ★, effet de latitude,
classement des codecs audité).

**Lancer les deux métriques et comparer les deux digests avant de toucher au
texte.** Si les verdicts diffèrent, c'est le résultat central du papier qui
change, pas un détail de présentation.

---

## Étape 6 — ce qu'il faut alors modifier dans `tex/tgrs_healpix_compress.tex`

Dans l'ordre du fichier :

1. **Abstract** — « seven real Sentinel-2 B04 scenes » → 40 régions ; remplacer
   tous les ratios et l'énoncé sur Holm ; retirer la phrase sur SZ3 comme
   meilleur codec général (voir l'audit ci-dessus) ; ajouter une phrase sur la
   référence de reconstruction commune si l'étape 4 est intégrée.
2. **Contributions** (3ᵉ puce) — « seven Sentinel-2 scenes » → 40 régions
   réparties en latitude ; ajouter l'effet de régime HEALPix comme contribution
   si `table_latitude_effect` donne un signal.
3. **§III-B** — préciser la taille de la fenêtre native et pourquoi : la boîte
   englobante de la tuile fait 901 px, mais l'ancre étant un point quelconque de
   la grille, la fenêtre garantie en est le double (1 801 px) — d'où 1 920 px
   halo compris. Et préciser que le patch 256² publié par
   healpix-resample en est le carré central.
4. **§IV-A Datasets** — remplacer `\input{tables/table_benchmark_scenes}` par
   `\input{tables/table_benchmark_regions}`, citer le DOI Zenodo
   `10.5281/zenodo.22210945` et le papier compagnon.
5. **§V** — remplacer le tableau γ★ écrit en dur (actuellement inline) par
   `\input{tables/table_gamma_star}`.
6. **§V-C** — supprimer entièrement le paragraphe « The `Wins (of 51)`
   denominator ... unresolved bookkeeping inconsistency » et réécrire le
   classement à partir de `table_codec_ranking` régénérée.
7. **Nouvelle sous-section Résultats** — effet de latitude / régime HEALPix,
   avec `table_latitude_effect` et `fig_latitude_effect`.
8. **Nouvelle sous-section Méthodologie** — la référence commune et le plancher
   de représentation (étape 4).
9. **Limitations** — retirer l'item « nine per-target significance tests are
   exploratory » si Holm passe désormais, retirer l'item « different
   within-representation reconstruction errors » si l'étape 4 est intégrée,
   retirer l'item « seven scenes, one tile each » ; remplacer l'item « oracle
   codec » par un renvoi à `table_fixed_codec` ; **ajouter** un item sur
   l'asymétrie de modèle Clough–Tocher / opérateur PSF.
10. **Conclusion** — réaligner sur le digest.

---

## Étape 7 — reproductibilité

Une fois les chiffres figés :

* transposer depuis healpix-resample `DATA_DEPOSIT.md`, `build_data_manifest.py`,
  `paper_data_guard.py`, `build_zenodo_archives.py` ;
* déposer les 40 fenêtres 1512² + les tuiles dérivées sur Zenodo, en renvoyant au
  DOI parent `10.5281/zenodo.22210945` ;
* ajouter `CITATION.cff`, `codemeta.json`, `.zenodo.json`, `git_commit.txt`.

---

## Ancres au bord d'un granule Sentinel-2

**Signature.** Tous les produits candidats, quelle que soit la date et quelle que
soit la tuile MGRS, sont écartés — et avant même le téléchargement, sur la seule
empreinte (« N rejected on footprint alone »). Un nuage varie d'une acquisition à
l'autre ; un rejet systématique signifie que la fenêtre demandée déborde du
polygone de données valides, qui est fixe pour une tuile MGRS et une orbite
données. Les granules Sentinel-2 ne dépassent leur carré MGRS nominal que de
quelques kilomètres, et ceux en bord de fauchée ne sont que partiellement
remplis.

**La vraie cause est le centrage, pas la taille.** Une fenêtre centrée sur
l'ancre doit faire le double de la boîte englobante de la tuile (1 801 px), parce
que l'ancre peut tomber n'importe où dans sa tuile, coin compris. Une fenêtre
centrée sur **la tuile qui contient l'ancre** n'a besoin que de la boîte
englobante plus le halo, soit environ 980 px — deux fois moins large. C'est
exactement la différence entre une fenêtre qu'un granule couvre et une fenêtre
qu'aucun granule ne couvre.

L'ancre continue de désigner le site et de sélectionner la tuile ; rien dans ce
calcul ne regarde une image, donc le choix ne peut favoriser aucune des deux
grilles.

**L'empreinte STAC classe, elle ne veto pas.** Le polygone de données valides est
une métadonnée grossière, parfois plus conservatrice que la donnée réelle. Les
produits dont l'empreinte couvre la fenêtre sont donc essayés en premier — cela
évite des téléchargements inutiles dans le cas normal — mais si aucun ne
convient, ceux que l'empreinte avait écartés sont lus quand même et c'est la
fraction de pixels réellement valides qui tranche. Le journal marque alors le
produit `+footprint-override` et l'affiche à l'écran. Laisser une métadonnée
seule éliminer une région serait la même erreur que de faire confiance
aveuglément à `items[0]`. Le patch 256² du papier compagnon n'est plus le carré central de la
fenêtre mais reste à l'intérieur de la tuile, autour de l'ancre — formulation à
reprendre telle quelle dans le papier.

**C'est désormais le comportement par défaut.** `build_region_windows.py` calcule
la fenêtre nécessaire région par région et la centre sur la tuile ; `--window` ne
sert plus qu'à forcer une taille, et `--centre anchor` à revenir à l'ancien
comportement (refusé sous 1 801 px).

```bash
python window_geometry.py                      # ce que chaque ancre exige
python build_region_windows.py --force         # les 40, fenêtre calculée par région
python build_tiles.py --force
```

Sur les 40 régions, la fenêtre centrée sur la tuile va de 984 à 1 194 px (médiane
1 135), contre 1 054 à 1 698 px centrée sur l'ancre — d'où le déblocage.

`window_geometry.py` compare aussi les tuiles déjà construites à celles que les
ancres sélectionnent. Les divergences attendues sont celles dont l'ancre est à
plus d'une demi-longueur de tuile (3 183 m) de son centre : le centre d'une tuile
voisine est alors plus proche de l'ancre que celui de sa propre tuile, et une
fenêtre centrée sur l'ancre retient la voisine. Les trois cas observés — iowa
(4,58 km), finland_boreal (4,67 km), lake_michigan (5,21 km) — sont exactement
les trois plus grands décalages du tableau, ce qui confirme que les deux calculs
sont cohérents.

Deux façons de finir, équivalentes scientifiquement :

* **uniforme** (recommandé, ~15 min de téléchargement) : ré-acquérir les 40 avec
  la règle unique ci-dessus, une seule phrase à écrire dans la section Données ;
* **minimal** : ne reprendre que les cinq régions concernées —
  `sao_paulo mato_grosso iowa finland_boreal lake_michigan` — en `--force` sur
  l'acquisition *et* sur le tuilage. Les 35 autres gardent exactement la même
  tuile ; seule la largeur de fenêtre diffère d'une région à l'autre, et elle
  n'est pas une grandeur mesurée.

**Si un site échoue encore** avec une fenêtre centrée sur sa tuile, il est
réellement à cheval sur un bord de granule permanent : la région sort de
l'échantillon. Critère de disponibilité de la donnée, évalué avant toute mesure
et identique pour les deux grilles, donc sans biais — mais à déclarer, avec le
décompte final, dans la section Données. Ne pas déplacer l'ancre pour sauver la
région : cela introduirait une sélection corrélée au site.

---

## Régions rejetées au tuilage : lire le diagnostic

Le message d'erreur donne maintenant, pour les quatre meilleures tuiles
candidates, le pourcentage d'inclusion dans la fenêtre et le pourcentage de
no-data dans l'empreinte. Trois cas :

**« none of the N candidate tiles fits ... best 98–99 % ».** C'était un défaut de
classement, corrigé : les candidates sont désormais rangées par centrage et toutes
essayées. Relancer `python build_tiles.py` suffit. Si le message persiste avec
toutes les candidates sous 100 %, la fenêtre est réellement trop petite pour ce
site (latitude élevée, forte rotation locale) : `--halo 16`, ou ré-acquérir avec
`--window 1800`.

**« no-data pixels ... above the 0.100 % ceiling » côté crop UTM.** Le même
plafond s'applique aux deux grilles. Sous le plafond, les quelques pixels isolés
du crop sont interpolés depuis leurs voisins finis — le pendant exact du retrait
des échantillons no-data de la triangulation HEALPix — et le compte est écrit
dans `tiles_index.csv` (`utm_nodata_filled`). Un détail à connaître : le lecteur
marque `DN == 0` comme no-data, ce qui sur de l'eau très sombre (lacs profonds,
ombres) signale parfois des pixels de réflectance réellement nulle plutôt que des
pixels manquants — d'où les 130 pixels du lac Taupo. Les comptes par région
restent donc auditables.

**« exceed the ... no-data ceiling » avec quelques centièmes de pour cent.**
Quelques pixels isolés. `--max-nodata 0.005` (0,5 %) reste défendable si vous le
déclarez dans le papier, mais vérifiez d'abord `nodata_in_tile_fraction` : le
plafond par défaut a été calibré pour laisser passer le bruit isolé et arrêter les
bords de granule.

**Fractions à deux chiffres (Amazon, Congo, São Paulo étaient à 21–32 %).** Ce
n'est pas du bruit : c'est un bord de granule ou un masque nuage dans la fenêtre.
Ne pas monter le plafond. Le bon geste est de changer de produit — l'acquisition
écarte maintenant d'elle-même tout produit dont la fenêtre a moins de
`--min-finite` (défaut 99,9 %) de pixels valides et passe à l'acquisition
suivante :

```bash
python build_region_windows.py --regions amazon congo sao_paulo --force
python build_tiles.py --regions amazon congo sao_paulo
```

Si aucune acquisition propre n'existe sur la fenêtre temporelle, la région sort de
l'échantillon — critère de disponibilité de la donnée, évalué avant toute mesure,
donc sans biais entre grilles, mais à déclarer dans la section Données.

---

## Les quatre régions sans COG public

**Symptôme** (rencontré au premier run) :

```
RasterioIOError: '/vsicurl/https://sentinel-s2-l2a.s3.us-west-2.amazonaws.com/
tiles/32/U/MU/2024/12/4/0/R10m/B04.jp2' does not exist ...
```

**Diagnostic.** Le bucket `sentinel-s2-l2a` est l'archive JP2 *requester-pays* de
Sinergise, pas le bucket public de COG. earth-search sert parfois l'asset `red`
depuis cette archive ; une lecture anonyme échoue alors, alors que l'item STAC est
parfaitement valide. Ce n'est ni transitoire ni propre à ce dépôt : la campagne
d'acquisition de `healpix-resample` a buté sur les mêmes quatre régions, avec les
mêmes hrefs, et c'est pour cela que ces quatre patches viennent du catalogue EOPF
dans le bundle publié.

**Ce que fait maintenant le script.** Il ne prend plus `items[0]` aveuglément. Pour
chaque région il construit une liste ordonnée de candidats — produit épinglé par le
manifest resample, puis même acquisition que le patch EOPF quand le bundle
resample est extrait localement (appariement sur tuile MGRS + jour, les
identifiants EOPF et earth-search ne se ressemblant pas), puis la recherche triée
par couverture nuageuse — sur trois collections (`sentinel-2-l2a`,
`sentinel-2-c1-l2a`, `sentinel-2-pre-c1-l2a`). Chaque candidat dont l'asset est un
JP2 ou vit dans l'archive requester-pays est écarté avant même la lecture, et une
`RasterioIOError` fait passer au suivant. Tout est journalisé.

**Relancer** — les 36 régions déjà acquises sont sautées, seules les 4 manquantes
sont retentées :

```bash
python build_region_windows.py
python build_region_windows.py --regions black_forest finland_boreal lake_constance lake_geneva --max-candidates 25
```

**Si elles échouent encore**, trois options, par ordre de préférence :

1. **Extraire le bundle Zenodo de healpix-resample** dans son checkout
   (`notebooks/data/multi_patch_sentinel2/`) puis relancer : le script lit alors
   `source_item_id` dans le store gelé et cherche la *même acquisition* sur
   earth-search. C'est la meilleure réponse en provenance — mêmes régions, mêmes
   acquisitions que le papier compagnon.
2. **Les fetcher depuis EOPF**, comme resample l'a fait, et fenêtrer à 1512 px
   soi-même. Même source exactement, mais code d'acquisition supplémentaire à
   écrire et convention radiométrique à revérifier (offset BOA).
3. **Continuer à n = 36.** C'est propre scientifiquement — le critère d'exclusion
   est la disponibilité de la donnée, évalué avant toute mesure, donc sans biais
   possible en faveur d'une grille — à condition de l'écrire noir sur blanc dans
   la section Données : 36 régions sur 40, les 4 écartées faute de COG public, et
   la différence de périmètre avec le papier compagnon assumée. Toute la chaîne
   d'analyse est pilotée par les données : `n = 36` passe partout sans
   modification.

Ne pas, en revanche, aller chercher ces 4 régions à une autre date ou une autre
position que les autres pour « sauver » le compte : cela introduirait une
sélection corrélée au site, ce que le protocole apparié ne pardonne pas.

---

## Dépannage

| Symptôme | Cause probable |
|---|---|
| `healpix-resample checkout not found` | poser `HEALPIX_RESAMPLE_ROOT` |
| `The AWS Access Key Id you provided does not exist` | credentials AWS dans l'environnement : `unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY` |
| `no level-10 tile is fully inside the window` | fenêtre trop petite ou ancre trop près d'un bord de granule : augmenter `--window`, ou accepter le rejet de la région |
| `the resampler kept N of 1048576 cells` | cellules hors enveloppe convexe : augmenter `--halo` |
| `--metric common needs the rmse_common column` | l'étape 4 n'a pas été lancée |
| `RasterioIOError ... B04.jp2 does not exist` | asset servi depuis l'archive requester-pays : voir « Les quatre régions sans COG public » |
| `none of the N candidate products yielded a readable public COG` | idem, après N tentatives : augmenter `--max-candidates`, ou accepter n = 36 |
| `none of the N candidate tiles fits` | voir « Régions rejetées au tuilage » |
| `granule footprint does not cover the ... window` / `rejected on footprint alone` | ancre au bord d'un granule : voir « Ancres au bord d'un granule Sentinel-2 » |
| `exceed the ...% no-data ceiling` | idem : bruit isolé → `--max-nodata`, bord de granule → ré-acquérir avec `--force` |
| `codec-ranking bookkeeping broken` | assertion volontaire : le décompte des victoires ne colle pas au nombre de points, ne pas publier la table |

---

## Ordre de grandeur du temps de calcul

À mesurer sur les deux premières régions, puis extrapoler : acquisition
~10-30 s/région, tuilage quelques minutes/région (Delaunay), balayage quelques
minutes/région pour les deux grilles, référence commune ~1-2 min/région en plus.
Compter une demi-journée machine pour la chaîne complète sur 40 régions, sans
surveillance grâce à la reprise automatique.
