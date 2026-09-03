# Brief : révision finale de `tex/tgrs_healpix_compress.tex`

Ce document sert de prompt/brief pour une nouvelle session dont l'unique objectif est de
corriger le papier IEEE GRSL (`tex/tgrs_healpix_compress.tex`) avec les résultats finaux,
validés, issus de l'investigation résumée ci-dessous. **Toutes les données sont maintenant
fiables et cohérentes entre elles** (notebook `test/stats_significance.ipynb` exécuté en
entier, "Restart & Run All", sans erreur, sur les CSV les plus récents). Le feu vert pour
éditer le papier est donné — contrairement aux échanges précédents où l'analyse devait
rester séparée de l'édition, cette fois le but explicite est de corriger le texte.

## 1. Ce qui a changé depuis la dernière version publiable du papier

Trois bugs méthodologiques majeurs ont été trouvés et corrigés, dans cet ordre :

1. **Overhead opérationnel UTM (γ=1.30) mal appliqué.** Une première tentative l'appliquait
   au mauvais terme (division du dénominateur brut), ce qui avantageait artificiellement
   UTM. La correction : γ multiplie directement les **octets compressés** du côté UTM (ou,
   de façon équivalente, le budget cible est divisé par γ avant interpolation quand on
   raisonne "RMSE à budget fixé"). C'est la convention utilisée partout maintenant.

2. **Découpage UTM non équivalent en aire.** L'ancien crop UTM était une bounding-box
   axis-aligned autour de l'empreinte en losange (tournée à 45°) de la tuile HEALPix, ce qui
   incluait ~2× trop de pixels (~900-940 px/côté au lieu de ~637). Corrigé par un carré UTM
   axis-aligned de **surface égale** (aire garantie identique car HEALPix est equal-area),
   côté = 1024 × Δeq(niveau 20) / 10m ≈ 637 px, centré sur le centroïde de la tuile HEALPix
   (plus proche pixel UTM). Voir `test/utm_codec_bench.ipynb` cellule 3 et la figure
   `fig_utm_crop_footprints.pdf`.

3. **Artefact de rééchantillonnage HEALPix (la découverte la plus importante).** Un spectre de
   puissance 2D montrait un excès de haute fréquence sur la tuile HEALPix, absent en UTM natif
   et absent des spectres tangent-plan du papier `healpix-resample`. Diagnostic complet :
   - Hypothèse écartée : déformation lente de la position des pixels (gradient lon/lat lisse,
     confirmé par diagnostic direct — pas de mécanisme haute fréquence possible).
   - Hypothèse intermédiaire (vraie mais secondaire) : les cellules HEALPix sont *equal-area
     mais pas equal-shape* — élongation locale mesurée jusqu'à ~1.7× entre les deux axes du
     raster aplati, ce qui réduit l'espacement effectif dans la direction "longue" à ~8.9 m,
     proche du pas natif Sentinel-2 (10 m) → risque de repliement de Nyquist réel, mais
     insuffisant à lui seul pour expliquer toute l'amplitude observée.
   - **Cause principale confirmée** : la génération des tuiles HEALPix utilisait
     `BilinearResampler` (`healpix_resample`), une interpolation de Shepard (IDW, Npt=4) dont
     le noyau `sigma_m` par défaut est calé sur l'échelle du pixel HEALPix (~6.2 m) et non sur
     le pas de la grille source UTM (~10 m) — un noyau trop étroit pour filtrer correctement,
     injectant du contenu haute fréquence parasite. **Remplacé par
     `scipy.interpolate.griddata(method='cubic')`** (interpolation Delaunay bicubique) dans
     `test/compress_test.ipynb` cellule 10. Confirmé par l'utilisateur : le spectre de
     puissance 2D redevient superposable point par point à celui d'UTM natif après ce
     changement — **c'était bien l'origine principale de l'écart de compressibilité observé**,
     pas une limite fondamentale de la géométrie HEALPix.

4. **Erreur d'agrégation entre scènes dans les Tables IV/V.** La première version calculait,
   *indépendamment à chaque cible RMSE/budget*, une médiane de tailles/RMSE *absolues* sur
   les 7 scènes, avec le "meilleur codec" re-choisi cible par cible — ce qui mélange des
   populations de scènes différentes selon la cible et casse la monotonicité (vérifié :
   `fig_rd_curves_per_scene.pdf`/`fig_utm_rd_curves_per_scene.pdf` par scène sont bien
   monotones individuellement). **Corrigé** : calcul du ratio HPX/UTM *à l'intérieur de
   chaque scène d'abord* (codec oracle par scène, cohérent avec le test de Wilcoxon apparié),
   puis médiane/range des 7 ratios. C'est désormais la méthodologie unique utilisée dans tout
   le notebook (Tables IV, V, et le test de significativité).

## 2. Résultats finaux validés (à utiliser tels quels dans le papier)

### Table III — ratio brut, sans compression (inchangée, non affectée par les bugs 3-4)
HPX raw / UTM opérationnel = **0.954** en moyenne (range 0.919–0.993), Wilcoxon p=0.016,
IC bootstrap 95% très serré autour de cette valeur. → HEALPix brut reste ~5% plus compact
qu'UTM opérationnel *avant toute compression*. Fichier source :
`test/results/storage_accounting.csv` / `tex/tables/table_storage_accounting.tex`.

### Table IV — ratio HPX/UTM à RMSE cible constante (par scène, médiane + [min,max])
| Cible RMSE | n | Ratio médian | Range | Favorise |
|---|---|---|---|---|
| 1e-4 | 7 | 1.59 | [1.26, 1.71] | UTM |
| 3e-4 | 7 | 1.47 | [0.61, 1.68] | UTM |
| 1e-3 | 7 | 1.24 | [0.94, 1.89] | UTM |
| 3e-3 | 5 | 0.80 | [0.54, 1.24] | HPX (non significatif, p=0.625) |
| 1e-2 | 5 | 0.91 | [0.68, 2.02] | HPX (non significatif, p=0.625) |

### Table V — ratio HPX/UTM de RMSE à budget de stockage constant (par scène)
| Budget (% partagé) | n | Ratio médian | Range | Favorise |
|---|---|---|---|---|
| 20% | 6 | 13.44 | [8.65, 33.51] | UTM (p=0.031) |
| 10% | 7 | 4.17 | [2.88, 4.43] | UTM (p=0.016) |
| 5% | 7 | 1.32 | [1.26, 1.44] | UTM (p=0.016) |
| 2% | 7 | 0.86 | [0.73, 0.96] | **HPX (p=0.016, significatif)** |

Source exacte : `test/results/significance_tests.csv`,
`tex/tables/table_storage_at_constant_rmse.tex`, `tex/tables/table_error_at_constant_size.tex`.

### Diagnostic clé : pourquoi le basculement de verdict (analyse "brut vs ajusté")
En retirant γ (ratio de compressibilité *intrinsèque*, sans overhead opérationnel), UTM reste
systématiquement plus efficace que HEALPix sur toute la plage testée, mais l'écart se réduit
progressivement avec la cible RMSE : **+107%, +91%, +61%, +4%, +18%** (de 1e-4 à 1e-2). Le
basculement de verdict dans les tables ci-dessus vient exactement de la comparaison entre cet
écart intrinsèque et les 30% de surcharge opérationnelle UTM : tant que l'écart brut dépasse
30% (les 3 cibles serrées), UTM reste devant même après retrait de sa propre surcharge ; dès
que l'écart brut tombe sous 30% (3e-3, 1e-2, et le budget 2%), l'absence de surcharge de
tuilage côté HEALPix devient l'argument décisif. **C'est un message central et défendable
pour le papier** : à haute fidélité, le raster HEALPix aplati reste intrinsèquement moins
compressible qu'un crop UTM natif (probablement lié à l'anisotropie de forme des pixels
equal-area) ; à compression agressive, cet écart intrinsèque devient plus petit que la
surcharge réelle de tuilage UTM, et HEALPix l'emporte *parce qu'il n'a pas cette surcharge*,
pas parce qu'il compresse mieux intrinsèquement.

### Table VI — classement des codecs (inchangée, non affectée)
SZ3 gagnant dans 45/51 points de fonctionnement HPX (rang moyen 1.47), JPEG 2000 lossy 9
victoires, JPEG XL lossy 5, DWT+zstd 2, ZFP 4. Toujours valide, aucune modification requise.

## 3. Figures à utiliser

- **`figures/fig_utm_crop_footprints`** (déjà inclus, `\label{fig:utm_crop_footprints}`,
  ligne ~136) : montre pour chaque scène le contour du carré UTM équivalent-aire (637×637)
  autour de l'empreinte HEALPix. **Vérifier qu'elle est bien référencée et expliquée dans le
  texte de la méthodologie (§III-B)** — c'est la figure qui rend concret tout le
  raisonnement du redécoupage équivalent-aire ; ne pas la laisser orpheline ou sous-utilisée.
- **`figures/fig_ratio_boxplots`** (nouvelle, générée dans `tex/figures/`,
  PAS ENCORE dans le papier) : box-plots des distributions de ratio HPX/UTM par scène, un
  panneau pour les cibles RMSE, un pour les budgets. **C'est le résumé visuel le plus direct
  du résultat central du papier** (transition UTM-favorisé → HPX-favorisé selon le régime de
  compression, avec la dispersion inter-scène visible). À ajouter comme figure principale des
  résultats, probablement en remplacement ou en complément de `fig_storage_at_constant_rmse`
  /`fig_error_at_constant_size` (l'ancienne figure 2 panneaux, ligne ~191, qui utilisait
  encore l'ancienne méthodologie de courbes agrégées en valeur absolue — à vérifier si elle
  reste pertinente ou doit être retirée/régénérée avec la nouvelle méthodologie par ratio).

## 4. Ce qu'il faut réécrire dans le papier, section par section

- **Abstract** (ligne ~54-56) : entièrement à réécrire. Le texte actuel affirme que HEALPix
  "does not inherently incur a storage penalty... across most of the practically relevant
  compression range" et que UTM ne devient meilleur "at the single most aggressive
  compression setting" — c'est l'**inverse** des résultats finaux. Nouvelle narrative : UTM
  reste plus efficace aux réglages de haute fidélité/fichier généreux (3 cibles RMSE serrées,
  3 budgets serrés), HEALPix devient compétitif voire meilleur aux réglages de compression
  agressive (2 cibles RMSE lâches — non significatif — et le budget le plus serré à 2% —
  significatif), et ce basculement s'explique par la comparaison entre l'écart de
  compressibilité intrinsèque (qui se réduit avec la cible) et la surcharge opérationnelle
  UTM de 30%.
- **Introduction, liste des contributions** (ligne ~73-80) : le 4e bullet ("HEALPix retains a
  directionally consistent... compression advantage... except at the single most aggressive
  compression setting") est à corriger dans le même sens. Ajouter une contribution sur la
  découverte méthodologique du resampler (sensibilité du résultat au noyau d'interpolation
  utilisé pour peupler la tuile HEALPix) — c'est un résultat scientifique en soi.
- **§III-A "From a Nested Chunk to a 2-D Raster"** (ligne ~122-127) : la phrase "This reshape
  is the only HEALPix-specific step in the pipeline; every compression step that follows
  operates on an ordinary 2-D float32 raster, exactly as it would on any other image" est à
  nuancer — on a montré que le choix du noyau de rééchantillonnage (et l'anisotropie de forme
  des pixels equal-area) a un effet de premier ordre sur la compressibilité, donc ce n'est pas
  aussi "ordinaire" qu'annoncé.
- **§III-B (méthodologie de crop UTM équivalent-aire)** : déjà rédigée correctement dans une
  passe précédente (637×637, dérivation exacte, référence à la figure) — vérifier qu'elle est
  cohérente avec la figure `fig_utm_crop_footprints` et ne nécessite pas d'autre changement
  de fond, seulement une relecture de cohérence globale.
- **Nouvelle sous-section méthodologique à ajouter** (dans §III ou en Discussion) : décrire
  brièvement l'épisode du resampler (BilinearResampler → griddata cubic), pourquoi ça compte,
  et ce que ça implique pour la reproductibilité (le choix du noyau d'interpolation utilisé
  pour peupler la grille HEALPix n'est pas un détail d'implémentation neutre).
- **§Results** : il y a un marqueur `\todo{}` juste après `\section{Results}` (ligne ~181,
  actuellement affiché en rouge dans le PDF) signalant que les tables/figures de cette
  section datent de l'ancienne méthodologie bounding-box. **Le retirer** et regénérer tout le
  texte de la section avec les chiffres de la partie 2 ci-dessus.
- **Discussion** (ligne ~211 et suivantes) : réécrire entièrement le paragraphe d'
  interprétation pour refléter la nouvelle histoire "régime-dépendant" plutôt que "HEALPix
  gagne presque partout". Intégrer le diagnostic brut-vs-ajusté (partie 2) comme explication
  mécanistique centrale.
- **Limitations** (ligne ~222) : ajouter un item sur la sensibilité au choix du noyau de
  rééchantillonnage (résultat pouvant varier significativement selon l'interpolateur utilisé
  pour peupler la grille HEALPix — démontré empiriquement) ; ajouter un item sur la taille
  d'échantillon réduite aux cibles les plus lâches (n=5/7, IQR large, résultat non
  significatif à ces cibles) ; garder/adapter l'item PSFresampler existant (le
  rééchantillonnage utilisé ici reste volontairement simple ; un futur pipeline PSF-aware
  ajoutera du détail fin réel et réduira l'efficacité de compression rapportée ici).
- **Conclusion** (ligne ~237) : réécriture complète pour correspondre exactement aux chiffres
  finaux et à la narrative régime-dépendante.

## 5. Fichiers de référence pour extraire les chiffres exacts

- `test/results/significance_tests.csv` — source de vérité pour tous les ratios/p-values.
- `test/results/storage_accounting.csv` — ratio brut (Table III).
- `tex/tables/table_storage_accounting.tex`, `table_storage_at_constant_rmse.tex`,
  `table_error_at_constant_size.tex`, `table_codec_ranking.tex` — déjà à jour, prêts à
  `\input`.
- `test/stats_significance.ipynb` — notebook canonique, seule source de vérité pour
  recalculer/vérifier n'importe quel chiffre.
- `tex/figures/fig_utm_crop_footprints.pdf`, `tex/figures/fig_ratio_boxplots.pdf` — figures
  finales prêtes à l'emploi.

## 6. Ce qui n'a PAS besoin d'être revérifié

Table VI (classement des codecs) n'a jamais dépendu des bugs 1-4 ci-dessus (comparaison
interne HPX uniquement, sans normalisation UTM) — laisser telle quelle. La bibliographie et
les citations ajoutées dans une passe précédente restent valides.
