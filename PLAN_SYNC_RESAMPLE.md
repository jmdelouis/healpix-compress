# Plan de synchronisation `healpix-compress` ↔ `healpix-resample`

État au 2 septembre 2026. Ce document remplace, pour la phase qui vient, la feuille de route
de `PAPER_REVISION_BRIEF.md` (qui reste valable pour l'édition purement rédactionnelle).

---

## 1. État du package `healpix-compress`

**Code (`hcc/`)** — 2 sous-paquets, ~20 modules, 13 fichiers de tests unitaires :

- `hcc/compression/` : implémentations « maison » (dct, decompose, grad_dct, grad_dwt,
  quadtree_wavelet, spiht, zerotree, entropy, zorder).
- `hcc/codecs/` : enrobage zarr v3 de ces algorithmes, exposés comme serializers via
  `project.entry-points."zarr.codecs"` (`hcc-dct`, `hcc-grad-dct`, `hcc-grad-dwt`,
  `hcc-decompose`).
- `hcc/compression/zorder.py` (`expand_tile`) est le seul module réellement utilisé par le
  papier : c'est le de-interleaving Morton de la Section III-A.

**Constat central : le package et l'article ont divorcé.** L'article annonce explicitement
« No custom or HEALPix-aware codec is used in this paper » et compare cinq codecs externes
(ZFP, SZ3, JPEG 2000, JPEG XL, DWT+zstd). Les codecs maison de `hcc` ne sont donc plus
l'objet du papier, et la couche zarr de `hcc` — qui est pourtant le livrable logiciel
intéressant (rendre ces codecs utilisables comme serializers zarr sur des chunks HEALPix
NESTED) — n'apparaît nulle part dans l'article. C'est le premier écart à résorber.

**Infrastructure** : `pyproject.toml` minimal (numpy/scipy/zarr/PyWavelets), pas de
dépendance à `healpix_resample` (les notebooks l'importent quand même), pixi.toml sans
lock, docs sphinx squelettiques (api.rst + index.rst), CI présente. Manquent, par rapport à
`healpix-resample` : `CITATION.cff`, `codemeta.json`, `.zenodo.json`, `data_manifest.csv`,
`DATA_DEPOSIT.md`, `DATA_LICENSE.md`, `ATTRIBUTION.md`, `REPRODUCING.md`, `pixi.lock`,
`git_commit.txt`.

**Notebooks (`test/`)** : `compress_test.ipynb` (acquisition + tuilage HEALPix),
`hpx_codec_bench.ipynb`, `utm_codec_bench.ipynb`, `stats_significance.ipynb` (notebook
canonique des chiffres du papier), plus `Untitled.ipynb` et `_adapters_extract.py` à
nettoyer. 22 CSV de résultats dans `test/results/`, 16 tables LaTeX, ~20 figures.

**Données actuelles** : 7 scènes Sentinel-2 L2A B04, fenêtres natives 1512×1512 px
(15,12 km), acquises à la volée via le STAC EOPF, stockées en zarr v2 non versionnées, sans
manifeste ni checksums. Variantes `_b08` présentes pour 8 stores (NDVI amorcé mais non
exploité). Aucun DOI, aucune provenance figée.

**Hygiène git** : ~30 fichiers non suivis (dont toutes les figures/tables finales),
3 répertoires `tex/_validate*` de compilation, `table_ndvi_stub_DO_NOT_PUBLISH.tex`.

---

## 2. État de l'article `tex/tgrs_healpix_compress.tex`

298 lignes, structure complète (Intro, Background, Méthodologie, Design expérimental,
Résultats, Discussion, Limitations, Conclusion, Annexe global-sphère, remerciements
GRID4EARTH). Le texte a déjà intégré la révision décrite dans `PAPER_REVISION_BRIEF.md` :
crop UTM équivalent-aire 637×637, overhead γ_tile=1.30 appliqué aux octets compressés,
agrégation par ratio intra-scène, correction de Holm, tableau γ★, section
III-C sur la sensibilité au noyau de rééchantillonnage.

**Ce qui est solide** : la formulation du problème, le protocole équivalent-aire, l'analyse
γ★ (écart intrinsèque vs overhead critique), la section Limitations (14 items,
inhabituellement honnête).

**Ce qui bloque une soumission** :

1. **Puissance statistique nulle.** n = 5 à 7 scènes, 9 tests, aucun ne survit à Holm
   (p ajusté >= 0,14). Le résultat central — le basculement de régime — est donc
   présente comme « exploratoire et directionnellement cohérent », ce qui est correct mais
   difficilement publiable comme résultat principal.
2. **Référence de reconstruction non commune** (Limitation 5) : la RMSE HEALPix est mesurée
   contre sa propre tuile déjà rééchantillonnée, la RMSE UTM contre le crop brut. Le papier
   lui-même désigne ce point, avec la sensibilité au noyau, comme le plus susceptible de
   changer les chiffres.
3. **Sensibilité au noyau non quantifiée** (Limitation 1) : l'épisode IDW → Clough–Tocher
   est raconté qualitativement, sans métrique de distance spectrale ni balayage de noyaux.
4. **Aucune reproductibilité des données** : pas de DOI, pas de manifeste, pas de commit
   figé. Rédhibitoire pour TGRS/GRSL en 2026.
5. **Incohérences résiduelles** : le dénominateur de `table_codec_ranking` (51) ne colle pas
   à la somme des victoires (65) — signalé dans le texte, ce qui est honnête mais ne peut
   pas rester tel quel ; `\markboth` dit TGRS alors que le brief parle de GRSL ; six entrées
   `.bib` non citées, dont `ieee_tgrs_scope` à supprimer impérativement.
6. **Bibliographie mince** (13 entrées citées) pour un full paper.

---

## 3. Ce que `healpix-resample` apporte

Le papier compagnon (`tex/main.tex`, 1408 lignes, quasi soumis) et son package fournissent
quatre choses directement réutilisables :

- **Un jeu de données figé et cité** : DOI Zenodo `10.5281/zenodo.22210945`, 4 scènes
  Sentinel-2 « cœur » + **40 patches régionaux** (4 classes x 10 régions), chacun avec
  latitude, longitude et identifiant produit S2 exact dans `data_manifest.csv`, plus
  SHA-256 par asset. Couverture en latitude réelle : Singapour (1,3 deg) → Finlande
  boréale, donc les deux régimes HEALPix (ceinture équatoriale |lat| < 41,8 deg et calottes
  polaires) sont peuplés.
- **Une famille de resamplers UTM → HEALPix** exposée proprement :
  `NearestResampler`, `BilinearResampler`, `BicubicResampler`, `CloughTocherResampler`,
  `KNeighborsResampler`, `PSFResampler`, plus `subset_for_parent_cell` et les utilitaires
  géométriques `fwhm_to_scale` / `cell_size_m` / `recommend_npt`.
- **Un opérateur direct** : `KNeighborsResampler.invert(hval)` (round-trip prédiction,
  Section III-E du papier resample) qui projette un champ HEALPix sur la grille source. Il
  n'existe pas d'équivalent dans `healpix-compress`.
- **Un modèle complet d'infrastructure de reproductibilité** : `DATA_DEPOSIT.md`,
  `build_data_manifest.py`, `paper_data_guard.py`, `build_zenodo_archives.py`,
  `publish_paper_assets.py`, `load_data_in_zenodo.ipynb`, `.zenodo.json`, `CITATION.cff`.
  Ces scripts sont paramétrables et transposables presque tels quels.

---

## 4. Le piège géométrique 256×256 (point signalé)

Les patches publiés par `healpix-resample` font **256×256 px natifs à 10 m, soit
2,56 × 2,56 km**. Une tuile porteuse HEALPix projetée en UTM est un **losange tourné**, donc
la fenêtre native nécessaire est le côté de la tuile multiplié par √2 (pire cas 45°).
Avec Δ_eq(20) = 6,2176 m :

| δ | cellules | côté tuile | carré UTM équi-aire | fenêtre native mini | taille brute f32 |
|------:|---------:|-----------:|--------------------:|--------------------:|-----------------:|
| 7  | 128 × 128   | 796 m   | 80 px  | 113 px | 64 KiB  |
| 8  | 256 × 256   | 1 592 m | 159 px | 226 px | 256 KiB |
| 9  | 512 × 512   | 3 183 m | 319 px | 451 px | 1 MiB   |
| 10 | 1024 × 1024 | 6 367 m | 637 px | **901 px** | 4 MiB |

*Correction apportée en cours de route :* la colonne « fenêtre native mini » ci-dessus
est la **boîte englobante de la tuile**, valable seulement si la fenêtre est centrée sur
la tuile. La fenêtre est en fait centrée sur l'ancre de la région, point quelconque de la
grille HEALPix, qui peut tomber sur une jonction à quatre tuiles : le centre de la tuile
la plus proche est alors à une demi-diagonale et son coin opposé à une demi-diagonale de
plus. La fenêtre garantie vaut donc le **double** de la boîte englobante — 1 801 px pour
δ = 10, d'où l'acquisition à 1 920 px halo compris.

Conclusion : **les patches 256 px ne peuvent pas porter la tuile δ = 10 du papier**. Le
plus grand δ admissible dans 256 px est δ = 8 (226 px de fenêtre requise, 15 px de
marge par bord).

Trois effets à ne pas subir passivement :

- **Effet de taille de tuile.** Passer de 4 MiB à 256 KiB de charge utile multiplie par 16 le
  poids relatif des en-têtes de codec. Au budget le plus serré (2 %), le budget passe de
  ~84 ko à ~5,2 ko, dont ~100 à 400 octets d'en-tête JPEG XL / JPEG 2000 : 2 à 8 % du
  budget, non négligeable et asymétrique entre codecs.
- **Effet de remplissage.** Si la tuile HEALPix déborde la fenêtre native, les cellules non
  couvertes sont remplies (NaN → moyenne, ou fill value). Ces zones se compressent
  quasi gratuitement et **gonflent artificiellement** le côté HEALPix. Le papier n'a
  actuellement qu'un affichage `Cover %` informatif ; il faut une assertion dure.
- **Effet de contenu effectif.** Une tuile « grande » mais peuplée à partir d'une source dont
  la résolution effective est plus grossière que le niveau porteur mesure la compressibilité
  de la redondance d'interpolation, pas celle de la donnée. C'est exactement la thèse du
  papier, donc le protocole doit garantir que le rapport entre pas natif et pas porteur reste
  fixe (10 m / 6,22 m = 1,608) pour toutes les scènes comparées.

**Décision retenue : ré-extraire des fenêtres 1512×1512 px aux 40 centres régions**, ce qui
préserve δ = 10, niveau 20 et le rapport 2,58× du papier, tout en héritant de la
provenance, de la couverture en latitude et de la conception d'échantillonnage de
`healpix-resample`.

---

## 5. Décisions retenues

| Sujet | Décision |
|---|---|
| Données | Ré-extraire des fenêtres 1512×1512 px aux 40 centres régions du manifest resample, mêmes produits S2 (lat/lon + `product=` dans `data_manifest.csv`). n : 7 → 40. |
| Géométrie | δ = 10 conservé (1024×1024 cellules niveau 20), crop UTM équivalent-aire 637×637. |
| Resampler | **Clough–Tocher conservé** comme opérateur de référence (pas de balayage de noyaux comme résultat ; la limitation reste ouverte, mais documentée). |
| Référence RMSE | **Référence commune UTM** : décodage HEALPix → `invert()` de `healpix_resample` → comparaison à la même observation UTM, avec plancher d'erreur de représentation reporté séparément. |

Conséquence à assumer explicitement : peupler la tuile en Clough–Tocher puis évaluer via
l'opérateur direct PSF de `healpix-resample` mélange deux modèles (la tuile n'est pas
l'inverse de l'opérateur utilisé pour la prédiction). Le plancher de représentation sera donc
plus élevé qu'avec `PSFResampler`, et il faudra le dire. Un run `PSFResampler` en simple
contrôle de sensibilité (même protocole, 40 régions) coûte peu et désamorcerait la question
d'un rapporteur ; il est proposé en option P3.4.

---

## 6. Plan par phases

### P0 — Hygiène préalable (0,5 jour)

1. Committer ou ignorer les ~30 fichiers non suivis ; supprimer `tex/_validate*`,
   `test/Untitled.ipynb`, `test/_adapters_extract.py`, `tex/*.log` du suivi.
2. Supprimer `table_ndvi_stub_DO_NOT_PUBLISH.tex` et toute référence.
3. Trancher la cible : TGRS (full paper, exige P4 + biblio étoffée) ou GRSL (letter, format
   plus court). Aligner `\markboth`, la classe et la longueur en conséquence.
4. Nettoyer le `.bib` : retirer `ieee_tgrs_scope` et les 5 autres entrées non citées, ou les
   citer réellement.

### P1 — Données communes (2 à 3 jours)

1. Écrire `test/build_region_windows.py` : lire `data_manifest.csv` de `healpix-resample`,
   pour chacune des 40 régions extraire une fenêtre 1512×1512 px B04 (et B08, cf. P4.3)
   centrée sur le même (lat, lon), depuis **le même produit** (`product=` du manifest) via
   earth-search ou EOPF selon la source enregistrée.
2. Garde-fous par région, tous journalisés dans le manifeste :
   - couverture valide de la fenêtre = 100 % (pas de no-data, nuage sous le seuil déclaré) ;
   - tuile parente niveau 10 dominante entièrement contenue dans la fenêtre
     (`Cover == 100,00 %`, assertion dure et non plus impression) ;
   - carré UTM 637×637 centré sur le centroïde de la tuile entièrement contenu ;
   - si échec, recentrage sur le candidat suivant ou élargissement, et trace de la raison.
3. Aligner les 4 scènes « cœur » (urban/water/forest/agriculture) : vérifier que les centres
   des scènes compress actuelles coïncident avec les patches resample ; sinon, adopter les
   centres resample et documenter le changement.
4. Décider du sort des 3 scènes actuelles absentes de resample (clouds, ocean, snow_ice) :
   soit les conserver comme scènes supplémentaires hors-manifest, soit les abandonner. Les
   conserver déforme la comparaison par classe ; recommandation : les garder en annexe
   uniquement.
5. Construire `data_manifest.csv` côté compress avec `build_data_manifest.py` adapté
   (asset_id, source, source_identifier, sha256, DOI parent, notebook consommateur).

### P2 — Pipeline de tuilage (1 à 2 jours)

1. Ajouter `healpix-resample` aux dépendances de `pyproject.toml`, avec une borne de version
   correspondant au commit/DOI du papier resample.
2. Extraire de `compress_test.ipynb` (cellule 10) un module `test/tiling.py` testable :
   choix de la tuile parente dominante, génération du bloc de 4^10 enfants,
   `CloughTocherResampler`, assertion de couverture, sauvegarde.
3. Remplacer `all_image_cell.npy` / `all_image_data.npy` (7 × 1M) par un store zarr indexé
   par région, avec attributs de provenance (région, classe, lat, lon, produit, resampler,
   version du package, commit).
4. Faire passer au moins un chemin de benchmark par la couche zarr de `hcc` (codec
   `hcc-*` comme serializer), afin que le package redevienne le livrable logiciel décrit
   dans l'article et non un dépôt parallèle.

### P3 — Référence de reconstruction commune (3 à 4 jours, poste le plus lourd)

1. Ajouter un module `test/common_reference.py` : pour chaque (région, codec, paramètre),
   décoder la tuile HEALPix, appliquer `KNeighborsResampler.invert()` pour prédire les
   observations sur la grille UTM native du crop 637×637, et calculer la RMSE contre la même
   observation UTM que celle utilisée côté UTM.
2. Reporter séparément le **plancher d'erreur de représentation** : RMSE de la tuile HEALPix
   non compressée, projetée par `invert()`, contre l'observation UTM. Toutes les courbes
   débit-distorsion HEALPix seront bornées par ce plancher ; c'est ce plancher qui rend la
   comparaison honnête.
3. Régénérer les tables IV/V et le tableau γ★ sous cette métrique commune, en
   conservant en annexe les valeurs sous métrique interne pour montrer l'écart entre les deux
   lectures.
4. *(Option recommandée)* rejouer le même protocole avec `PSFResampler` (configuration du
   papier resample : FWHM 12,5 m, s = 10,62 m, q = 22, λ = 1e-3) comme unique contrôle
   de sensibilité au noyau, sans en faire un balayage complet. Coût marginal une fois P3.1
   écrit ; ferme partiellement la Limitation 1.

### P4 — Analyses rendues possibles par n = 40 (2 à 3 jours)

1. **Puissance statistique** : rejouer `stats_significance.ipynb` avec 40 paires. Wilcoxon
   exact apparié sur 40 échantillons donne une résolution très supérieure ; la correction de
   Holm sur 9 cibles redevient franchissable. Reporter des intervalles bootstrap plutôt que
   de simples médianes.
2. **Effet de latitude (résultat nouveau).** Les cellules HEALPix sont equal-area mais pas
   equal-shape, et l'anisotropie locale diffère entre ceinture équatoriale (|lat| < 41,8°)
   et calottes. Le papier invoque déjà une élongation jusqu'à ~1,7× sans la relier à la
   position. Avec 40 régions réparties en latitude, tester la corrélation entre
   |latitude| (ou le régime équatorial/calotte) et l'écart intrinsèque HPX/UTM. C'est
   l'apport scientifique le plus direct de la synchronisation, et il n'exige aucune donnée
   supplémentaire au-delà de P1.
3. **NDVI** : si B08 est extrait en P1, calculer l'effet de la compression sur le NDVI aux
   points de fonctionnement appariés, et remplacer l'argument indicatif actuel des
   Limitations par une mesure. C'est la réserve n° 1 identifiée dans
   `review_soumission_tgrs.md` pour une soumission TGRS.
4. **Analyse par classe** : 10 régions par classe permettent de reporter la comparaison
   par type de couverture avec une dispersion, et non plus une scène unique par classe.
5. Auditer et corriger le dénominateur de `table_codec_ranking` (51 vs 65) contre le
   notebook source ; vérifier que SZ3 n'est pas simplement absent des cibles les plus
   agressives par manque de couverture de balayage.

### P5 — Rédaction (2 à 3 jours)

1. Section Données : remplacer la description « sept scènes » par le protocole 40 régions,
   citer le DOI Zenodo de `healpix-resample` et le papier compagnon, décrire la fenêtre
   1512 px et le critère de couverture 100 %.
2. Nouvelle sous-section méthodologique : la référence de reconstruction commune via
   l'opérateur direct, avec le plancher de représentation.
3. Section Résultats : réécrire avec les chiffres n = 40 ; le paragraphe sur la correction de
   Holm change de nature (de « rien ne survit » à un énoncé de significativité réelle, dans
   un sens qu'il faudra accepter tel qu'il sortira).
4. Ajouter la sous-section « effet de latitude » si P4.2 donne un signal.
5. Limitations : retirer les items 4 (puissance), 5 (référence commune), et atténuer l'item 1
   si P3.4 est réalisé ; ajouter un item sur la cohérence de modèle (Clough–Tocher peuplant la
   tuile, opérateur PSF pour la prédiction).
6. Positionnement explicite des deux papiers : `healpix-resample` définit l'opérateur de
   transfert vers HEALPix, `healpix-compress` en mesure le coût de stockage. Citation croisée
   dans les deux sens (côté resample, la section « Implications for HEALPix EO Data Cubes »
   évoque déjà « storage and compression design » et peut porter la référence avant).
7. Étoffer la bibliographie : compression de cartes HEALPix en astronomie (CMB/Planck),
   formats cloud-native (COG, Zarr, STAC), comparatifs de compression d'imagerie satellite.

### P6 — Reproductibilité et dépôt (1 à 2 jours)

1. Transposer depuis `healpix-resample` : `DATA_DEPOSIT.md`, `build_data_manifest.py`,
   `paper_data_guard.py`, `build_zenodo_archives.py`, `REPRODUCING.md`.
2. Ajouter `CITATION.cff`, `codemeta.json`, `.zenodo.json`, `DATA_LICENSE.md`,
   `ATTRIBUTION.md`, `git_commit.txt`, `pixi.lock`.
3. Publier un dépôt Zenodo compagnon : fenêtres 1512 px aux 40 centres + tuiles HEALPix
   dérivées, avec renvoi explicite au DOI parent `10.5281/zenodo.22210945` et aux mêmes
   identifiants produits Copernicus (vérifier les termes de redistribution Sentinel-2, déjà
   traités côté resample dans `DATA_LICENSE.md`).
4. Mode `OFFLINE=True` par défaut dans les notebooks, comme côté resample.

---

## 7. Ordre et charge

| Phase | Charge | Depend de |
|---|---|---|
| P0 hygiène | 0,5 j | — |
| P1 données | 2 à 3 j | P0 |
| P2 pipeline | 1 à 2 j | P1 |
| P3 référence commune | 3 à 4 j | P2 |
| P4 analyses n=40 | 2 à 3 j | P3 |
| P5 rédaction | 2 à 3 j | P4 |
| P6 dépôt | 1 à 2 j | P1, en parallèle de P5 |

Total : environ 12 à 17 jours de travail effectif. P3 est le poste critique et le seul qui
puisse remettre en cause les chiffres publiés ; il doit etre fait avant toute reecriture de
la section Résultats.

---

## 8. Points de vigilance

- **Ne pas rejouer le piège de couverture.** Toute tuile dont la couverture native n'est pas
  à 100 % doit être rejetée, pas remplie. Un remplissage silencieux inverserait le résultat
  en faveur de HEALPix pour de mauvaises raisons.
- **Le rapport 10 m / 6,22 m doit rester constant** sur toutes les régions : c'est la
  grandeur dont le papier mesure le coût. Toute région où le produit S2 n'est pas à 10 m
  natif doit être écartée.
- **La comparaison par codec oracle reste optimiste** (Limitation 6). Avec 40 régions, la
  comparaison à codec fixé devient statistiquement lisible : la produire en même temps.
- **Le changement de référence RMSE peut inverser des conclusions.** Prévoir de publier les
  deux lectures et de ne pas figer le récit avant P3.3.
- **Ne pas étendre le périmètre** : multi-bandes, multi-capteurs, sphère globale et
  balayage complet de noyaux restent hors de cette passe.
