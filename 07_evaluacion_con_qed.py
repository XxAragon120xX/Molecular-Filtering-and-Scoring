import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen, QED
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from io import BytesIO

# Función para calcular propiedades moleculares incluyendo QED
# Función para calcular propiedades moleculares incluyendo QED
def calcular_propiedades(mol):
    if mol is not None:
        qed = QED.default(mol)
        qed_props = QED.properties(mol)  # Obtener las propiedades de QED
        return {
            'TPSA': Descriptors.TPSA(mol),
            'MW': Descriptors.ExactMolWt(mol),
            'nHA': Descriptors.NumHAcceptors(mol),
            'nHD': Descriptors.NumHDonors(mol),
            'LogP': Crippen.MolLogP(mol),
            'SA': Descriptors.TPSA(mol) / Descriptors.ExactMolWt(mol),  # Aproximación simple
            'QED': qed,
            'QED_MW': qed_props.MW,      # Peso molecular en las propiedades de QED
            'QED_ALOGP': qed_props.ALOGP,
            'QED_HBA': qed_props.HBA,
            'QED_HBD': qed_props.HBD,
            'QED_PSA': qed_props.PSA,
            'QED_ROTB': qed_props.ROTB,
            'QED_AROM': qed_props.AROM,
            'QED_ALERTS': qed_props.ALERTS
        }
    return {}


# Función para evaluar la idoneidad de una molécula
def evaluar_molecula(props, standard_value):
    score = 0
    # Propiedades que queremos bajas
    for prop in ['TPSA', 'MW', 'nHA', 'nHD', 'SA']:
        score += 1 / (props[prop] + 1)  # +1 para evitar división por cero
    
    # QED (más alto es mejor)
    score += props['QED'] * 10  # Damos más peso al QED
    
    # Standard Value (más bajo es mejor)
    score += 1 / (standard_value + 1)
    
    return score

# Cargar el archivo CSV
archivo_csv = "archivo_filtrado.csv"  # Cambia esto al nombre de tu archivo CSV
df = pd.read_csv(archivo_csv)

# Convertir la columna "Standard Value" a tipo numérico, ignorando los errores
df["Standard Value"] = pd.to_numeric(df["Standard Value"], errors='coerce')

# Eliminar filas con valores de "Standard Value" NaN
df = df.dropna(subset=["Standard Value"])

# Calcular propiedades y evaluar cada molécula
mejor_molecula = None
mejor_score = -float('inf')
mejor_props = None
mejor_standard_value = None
mejor_smiles = None

for index, fila in df.iterrows():
    smiles = fila["Smiles"]
    standard_value = fila["Standard Value"]
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        props = calcular_propiedades(mol)
        score = evaluar_molecula(props, standard_value)
        if score > mejor_score:
            mejor_molecula = mol
            mejor_score = score
            mejor_props = props
            mejor_standard_value = standard_value
            mejor_smiles = smiles

# Visualizar la mejor molécula
if mejor_molecula:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    
    # Dibujar la molécula
    img = Draw.MolToImage(mejor_molecula, size=(300, 300))
    img_byte_arr = BytesIO()
    img.save(img_byte_arr, format='PNG')
    img_byte_arr.seek(0)
    img = mpimg.imread(img_byte_arr)
    
    ax1.imshow(img)
    ax1.axis('off')
    ax1.set_title("Mejor Molécula")

    # Mostrar las propiedades
    props_text = "Propiedades:\n\n"
    for prop, value in mejor_props.items():
        props_text += f"{prop}: {value:.4f}\n"
    props_text += f"\nStandard Value: {mejor_standard_value:.4f}"
    props_text += f"\nSMILES: {mejor_smiles}"
    
    ax2.text(0.1, 0.5, props_text, fontsize=10, verticalalignment='center')
    ax2.axis('off')

    plt.tight_layout()
    plt.show()
else:
    print("No se encontró ninguna molécula válida.")

print(f"Mejor SMILES: {mejor_smiles}")
print(f"Mejor Standard Value: {mejor_standard_value}")
print(f"Mejor QED: {mejor_props['QED']}")