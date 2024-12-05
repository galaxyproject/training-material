---
title:  Importer via un lien
area: data upload
box_type: tip
layout: faq
contributors: [yvanlebras]
---

* Copier le lien
* Ouvrez le gestionnaire de téléchargement Galaxy ({% icon galaxy-upload %} en haut à droite du panneau d'outils)
{% if include.reset_form %}
* Selectionnez le bouton **Reset** en bas de ce formulaire. Si le bouton est grisé, passez à l'étape suivante.
{% endif %}
{% if include.collection %}
* Selectionnez **Collection** ci-dessous
{% endif %}
{% if include.collection_type %}
* Selectionnez **Type de collection** et choisir `{{ include.collection_type }}`
{% endif %}
* Selectionnez **Coller/Récupérer les données**
* Collez le lien dans le champ de texte
{% if include.link %}
  `{{ include.link }}`
{% endif %}
{% if include.link2 %}
  `{{ include.link2 }}`
{% endif %}
{% if include.format %}
* Modifier **Type (tout définir) :** de "Détection automatique" à `{{ include.format }}`
{% endif %}
{% if include.genome %}
* Remplacez **Genome** par `{{ include.genome }}`
{% endif %}
* Appuyez sur Start**
{% if include.collection %}
* Cliquez sur **Construire** lorsqu'il est disponible
{% if include.pairswaptext %}
* Assurez-vous que les "reads" avant et arrière sont définies sur {{ include.pairswaptext }}, respectivement.
    * Cliquez sur **Swap** sinon
{% endif %}
* Une convention de dénomination utile consiste à utiliser
{% if include.collection_name_convention %}
    * Es útil usar la próxima manera de nombrar archivos: {{ include.collection_name_convention }}
{% endif %}
{% if include.collection_name %}
    * {{ include.collection_name }}
{% endif %}
* Cliquez sur **Créer une liste** (et attendez un peu)
{% else %}
* Ferme la fenêtre
{% endif %}
{% if include.renaming == undefined or include.renaming == true %}
* Galaxy utilise les URL comme noms par défaut, vous devrez donc les remplacer par des URL plus utiles ou informatives.
the window
{% endif %}
