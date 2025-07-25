import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from io import BytesIO
import random

# Función para calcular propiedades moleculares
def calcular_propiedades(mol):
    if mol is not None:
        return {
            'TPSA': Descriptors.TPSA(mol),
            'MW': Descriptors.ExactMolWt(mol),
            'nHA': Descriptors.NumHAcceptors(mol),
            'nHD': Descriptors.NumHDonors(mol),
            'LogP': Crippen.MolLogP(mol),
            'SA': Descriptors.TPSA(mol) / Descriptors.ExactMolWt(mol),  # Aproximación simple
            'MDCK': random.uniform(0, 100),  # Simulado
            'BBB': random.uniform(0, 1),  # Simulado
            'F': random.uniform(0, 100),  # Simulado
            'CYP2C9_T12': random.uniform(0, 24),  # Simulado
            'hERG': random.uniform(0, 100),  # Simulado
            'ROA': random.choice(['Oral', 'Intravenoso', 'Tópico'])  # Simulado
        }
    return {}

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
fig, axs = plt.subplots(1, 3, figsize=(20, 10))

# Iterar sobre las muestras y mostrar cada molécula
for ax, (index, fila) in zip(axs, muestras.iterrows()):
    # Obtener la estructura SMILES, el valor estándar y el Molecule ChEMBL ID
    smiles = fila["Smiles"]
    valor = fila["Standard Value"]
    ch_emb_id = fila["Molecule ChEMBL ID"]
    
    # Convertir SMILES a una molécula RDKit
    mol = Chem.MolFromSmiles(smiles)
    
    # Calcular propiedades
    props = calcular_propiedades(mol)
    
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
        
        # Crear el título con las propiedades
        title = f"ID: {ch_emb_id}\nValue: {valor:.2f}\n"
        for prop, value in props.items():
            if isinstance(value, float):
                title += f"{prop}: {value:.2f}\n"
            else:
                title += f"{prop}: {value}\n"
        
        ax.set_title(title, fontsize=8)

# Ajustar el espacio entre subgráficas
plt.tight_layout()
plt.show()