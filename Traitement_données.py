#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector, Button

# --- 1. Importation des données ---
fichier = "DATA/adrian0.txt"

df = pd.read_csv(fichier, 
                 sep="\t", 
                 skiprows=12, 
                 engine="python")

df.columns = ["abs time(s)", "Fz(N)"]

# --- 2. Création du graphique ---
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(df["abs time(s)"], df["Fz(N)"], color="blue", label="Fz(N)")
ax.set_xlabel("Temps (s)")
ax.set_ylabel("Fz (N)")
ax.set_title("Sélectionne la zone à garder avec la souris")
ax.legend()

# --- 3. Variables globales pour stocker la sélection ---
selection = {}

# --- 4. Fonction appelée lors de la sélection ---
def onselect(xmin, xmax):
    selection["xmin"] = xmin
    selection["xmax"] = xmax
    print(f"Sélection : de {xmin:.3f}s à {xmax:.3f}s")
    ax.axvspan(xmin, xmax, color='gray', alpha=0.3)
    fig.canvas.draw_idle()

# --- 5. Création du SpanSelector ---
span = SpanSelector(ax, onselect, 'horizontal', useblit=True,
                    props=dict(facecolor='gray', alpha=0.3))

# --- 6. Fonction pour sauvegarder la sélection ---
def save_selection(event):
    if "xmin" in selection and "xmax" in selection:
        xmin, xmax = selection["xmin"], selection["xmax"]
        df_trimmed = df[(df["abs time(s)"] >= xmin) & (df["abs time(s)"] <= xmax)].copy()
        df_trimmed["abs time(s)"] -= df_trimmed["abs time(s)"].iloc[0]  # recaler le temps à 0

        # Affichage du résultat
        print(f"{len(df_trimmed)} lignes conservées.")
        print(df_trimmed.head())

        # Sauvegarde dans un nouveau fichier
        output_file = "DATA_correction/adrian0_trimmed.csv"
        df_trimmed.to_csv(output_file, index=False)
        print(f"✅ Données sauvegardées dans {output_file}")

        # Fermer la figure
        plt.close(fig)
    else:
        print("⚠️ Aucune sélection effectuée.")

# --- 7. Bouton Sauvegarder ---
ax_save = plt.axes([0.8, 0.05, 0.1, 0.075])
btn_save = Button(ax_save, "Sauvegarder")
btn_save.on_clicked(save_selection)

plt.show()
