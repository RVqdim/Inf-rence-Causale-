# Inf-rence-Causale-

Voici le repository qui rassemble les travaux réalisés pendant le projet :


## Simulations : 

Le fichier "simulations" contient l'application streamlit dans **causalite_simu_sensibility.py** et le nettoyage des données dans **data_cleaning.py**

Le dossier ***data*** reçois le fichier **df_filtered.csv** en sortie de **data_cleaning.py** 

Par soucis de stockage, les données ne sont pas présentes dans le dossier data. Pour les récupérer, il faut d'abord mettre le fichier genotype complet dans ***data***, le renommer **genotype.csv**, puis faire tourner **data_cleaning.py**. Le fichier **df_filtered** sera créé et il sera possible de faire tourner l'app.

Pour faire tourner l'app, ouvrir un terminal et écrire : *streamlit run causalite_simu_sensibility.py*


