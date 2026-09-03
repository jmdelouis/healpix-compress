# Revue de fond — "Compressing the HEALPix Resolution Mismatch for Earth Observation Data"

## Verdict

Le papier est solide, honnête intellectuellement, et l'observation centrale (95,4 %, 91,9–99,3 % sur les 7 scènes, systématiquement <100 % avant toute compression) est un résultat réel et bien établi. Mais en l'état, ce n'est pas encore un article TGRS "full paper" soumissible sans risque : la portée expérimentale (1 capteur, 1 bande, 7 scènes, 1 tuile par scène) et une formulation du résultat de compression ("statistiquement indiscernable") plus forte que ce que les données démontrent réellement sont les deux points qu'un reviewer va relever en premier. Ce sont des corrections de quelques jours de travail, pas une refonte — le papier n'a pas besoin d'une nouvelle campagne de calcul pour être défendable, juste d'un ajustement des affirmations et de deux ajouts (test statistique, related work).

## Ce qui est solide

Le renversement conceptuel — le mismatch de résolution HEALPix n'est pas une perte d'information mais de la redondance compressible — est formulé avec précision, sans sur-affirmation ("we deliberately avoid claiming... band-limited"), ce qui est exactement la rigueur attendue en review. La méthodologie de crop UTM "footprint-exact" est une vraie contribution méthodologique, correctement mise en défaut dans les limitations (biais conservateur du bounding box, qui joue *contre* la thèse du papier — donc crédible). La section Limitations est inhabituellement complète et honnête pour un draft d'auteur (7 points, dont plusieurs qui affaiblissent directement le résultat plutôt que de l'enjoliver) : c'est un vrai atout en review, les rapporteurs le remarquent favorablement. Le choix d'abandonner le codec maison et de l'assumer explicitement en Discussion ("closed part of the gap... but never matched them") est également un point de crédibilité plutôt qu'une faiblesse.

## Ce qui doit être corrigé avant soumission

**"Statistiquement indiscernable" n'est pas démontré.** L'abstract et la liste de contributions emploient ce terme, mais aucun test statistique n'est réalisé — seulement une lecture qualitative de "quelques points de pourcentage d'écart, l'avantage change de camp". Avec n=7 scènes, un test apparié (Wilcoxon signé, ou au moins un intervalle de confiance bootstrap sur l'écart médian) est faisable et rendrait l'affirmation défendable telle quelle. Sinon, il faut reformuler partout en "empiriquement comparables" / "dans une marge de quelques points de pourcentage" — un reviewer statisticien retoquera "statistically indistinguishable" sans test à l'appui, quasi systématiquement.

**Revue de littérature trop mince pour un full paper TGRS.** La bibliographie utilisée compte 13 entrées citées (5 codecs, HEALPix, H3, S2, OGC/ISO DGGS). Pour TGRS, une fourchette de 30 à 50 références est typique. Manquent notamment : des travaux antérieurs sur la compression de cartes HEALPix en astronomie (le domaine d'origine de HEALPix a une littérature dédiée sur la compression de cartes du CMB/Planck) ; des références sur les formats cloud-native géospatiaux (COG, Zarr, STAC) puisque l'introduction s'appuie explicitement sur cet argument ; des travaux comparatifs de compression d'imagerie satellite existants (au-delà des papiers fondateurs de chaque codec). Le fichier `.bib` contient d'ailleurs 6 entrées non citées dans le texte (`esa_sentinel2_msi`, `zarr_v3_spec`, `gdal_zarr`, `healpix_geo_docs`, `papadopoulos2016tiledb`, et surtout `ieee_tgrs_scope` — une citation de la page "Aim and Scope" de la revue elle-même, à supprimer, cela ne doit jamais apparaître dans une bibliographie soumise) : soit à citer réellement (contexte cloud-native, format Zarr), soit à retirer du `.bib`.

**La boucle motivation → résultats n'est pas refermée.** L'introduction justifie DGGS/HEALPix en partie par l'accès local efficace et le chunking pour le cloud ("efficient local access", "chunk contiguity controls I/O performance"). Mais les résultats ne reviennent jamais sur l'accès partiel/range-read : aucun des 5 codecs n'est comparé sous cet angle alors que c'est un argument moteur du papier. Soit retirer cet argument de l'intro, soit ajouter une phrase de discussion dessus (même qualitative : JPEG 2000/JPEG XL supportent le décodage progressif/tuilé, SZ3/ZFP moins nativement).

**Portée expérimentale à assumer plus frontalement dans le titre/abstract, ou à élargir.** Une seule bande (B04), un seul capteur (Sentinel-2), une tuile par scène. C'est déjà listé en Limitations — bien — mais un rapporteur TGRS demandera probablement soit un second capteur/bande avant acceptation, soit au minimum que l'abstract signale explicitement "single-band, single-sensor study" pour ne pas créer d'attente déçue à la lecture du titre général ("Earth Observation Data").

**Absence de produit dérivé (NDVI, etc.).** Déjà noté en Limitations comme travail futur — correct de le signaler, mais pour TGRS (orienté applications), l'absence totale de lien avec un produit scientifique réel est probablement la plus grosse réserve d'un rapporteur, plus que tout point méthodologique. À minima, une phrase de discussion sur l'impact attendu (ou non) de la RMSE observée sur un indice comme NDVI aiderait, même sans nouvelle expérience.

## Point méthodologique à vérifier (pas rédactionnel)

Le Tableau `codec_ranking` classe les codecs sur des points où chacun n'est évalué que dans sa plage de RMSE atteignable (pas d'extrapolation, cf. §Protocol). Il faut s'assurer que SZ3 n'est pas simplement absent des cibles les plus agressives (RMSE=1e-2, budget=2%) faute d'avoir un point de sweep qui va jusque-là — auquel cas JPEG 2000/JPEG XL "gagneraient" par forfait plutôt que par supériorité réelle à ces cibles. Si c'est le cas, il faut le dire explicitement (ce serait cohérent avec l'explication structurelle donnée, mais la nuance "gagne parce que seul évaluable" vs "gagne parce que meilleur" change la force de la conclusion en §V-C).

## Verdict par niveau de venue

- **IEEE TGRS (full paper), en l'état** : risque réel de "major revision" pour portée expérimentale limitée + affirmation statistique non étayée + related work insuffisant. Pas rédhibitoire, mais probable.
- **IEEE GRSL (Letters)** : la portée actuelle (une question précise, une réponse claire, méthodologie propre, résultat honnête) correspond bien au format court. À considérer sérieusement si l'objectif est une publication rapide et à moindre risque de rejet.
- **TGRS après corrections listées ci-dessus** (test stat + related work + phrase sur random-access + reformulation portée) : redevient un candidat raisonnable, le fond scientifique le supporte.

## Checklist avant soumission

1. Ajouter un test statistique (Wilcoxon apparié ou IC bootstrap) pour la comparaison post-compression, ou reformuler "statistically indistinguishable" → "empirically comparable, within a few percentage points" partout (abstract, contributions, conclusion).
2. Étoffer la bibliographie (cloud-native geospatial, compression HEALPix/CMB en astronomie) et retirer les 6 entrées `.bib` non citées, en particulier `ieee_tgrs_scope`.
3. Ajouter une phrase de discussion sur l'accès partiel/range-read pour reboucler avec l'argument cloud-native de l'introduction.
4. Vérifier si SZ3 est absent des cibles les plus agressives par manque de couverture de sweep plutôt que par infériorité réelle, et le préciser si c'est le cas.
5. Expliciter "single-band, single-sensor" dans l'abstract, pas seulement en Limitations.
6. Optionnel mais utile : une phrase (même sans nouvelle expérience) sur l'impact attendu de la RMSE mesurée sur un produit dérivé comme le NDVI.
