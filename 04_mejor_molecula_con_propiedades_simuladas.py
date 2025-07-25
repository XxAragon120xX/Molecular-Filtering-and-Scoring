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

# Seleccionar la molécula con el valor más bajo de "Standard Value"
fila = df.loc[df["Standard Value"].idxmin()]

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
    plt.figure(figsize=(8, 6))
    plt.imshow(img)
    plt.axis('off')
    
    # Mostrar las propiedades al lado de la imagen
    props_str = "\n".join([f"{prop}: {value:.2f}" if isinstance(value, float) else f"{prop}: {value}" for prop, value in props.items()])
    plt.text(1.05, 0.5, props_str, fontsize=12, va='center', transform=plt.gca().transAxes)
    
    # Título con el ID y el valor estándar
    plt.title(f"ID: {ch_emb_id}\nValue: {valor:.2f}", fontsize=14)
    
    plt.show()
