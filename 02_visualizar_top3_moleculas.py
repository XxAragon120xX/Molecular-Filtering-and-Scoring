import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from io import BytesIO

# Cargar el archivo CSV filtrado
archivo_csv = "archivo_filtrado.csv"  # Cambia esto al nombre de tu archivo CSV

# Leer el archivo CSV
df = pd.read_csv(archivo_csv)

# Convertir la columna "Standard Value" a tipo numérico, ignorando los errores
df["Standard Value"] = pd.to_numeric(df["Standard Value"], errors='coerce')

# Eliminar filas con valores de "Standard Value" NaN
df = df.dropna(subset=["Standard Value"])

# Seleccionar las 3 filas con los valores más bajos de "Standard Value"
muestras = df.nsmallest(3, "Standard Value")

# Crear una figura y un conjunto de ejes para mostrar las moléculas
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Iterar sobre las muestras y mostrar cada molécula
for ax, (index, fila) in zip(axs, muestras.iterrows()):
    # Obtener la estructura SMILES, el valor estándar y el Molecule ChEMBL ID
    smiles = fila["Smiles"]
    valor = fila["Standard Value"]
    ch_emb_id = fila["Molecule ChEMBL ID"]
    
    # Convertir SMILES a una molécula RDKit
    mol = Chem.MolFromSmiles(smiles)
    
    # Dibujar la molécula y convertir la imagen a un formato adecuado para matplotlib
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        img_byte_arr = BytesIO()
        img.save(img_byte_arr, format='PNG')
        img_byte_arr.seek(0)
        img = mpimg.imread(img_byte_arr)
        
        # Mostrar la imagen de la molécula
        ax.imshow(img)
        ax.axis('off')
        ax.set_title(f"ID: {ch_emb_id}\nValue: {valor}")

# Ajustar el espacio entre subgráficas
plt.tight_layout()
plt.show()

