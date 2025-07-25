import pandas as pd

# Cargar el archivo CSV original
archivo_csv = "Base de datos.csv"  # Cambia la ruta al archivo CSV original

# Leer el archivo CSV con pandas
df = pd.read_csv(archivo_csv, delimiter=';')

# Seleccionamos solo las columnas "Smiles", "Standard Value", "Standard Units"
df_filtrado = df[["Molecule ChEMBL ID","Smiles", "Standard Value", "Standard Units"]]

# Guardar el nuevo archivo CSV con las columnas seleccionadas
nuevo_csv = "moleculas filtrada.csv"
df_filtrado.to_csv(nuevo_csv, index=False)

print(f"Nuevo archivo CSV generado: {nuevo_csv}")

#////////////////////////////////////////////////////////////////////////////////////////////////////



# Cargar el archivo CSV filtrado previamente
segudno_archivo_csv = "moleculas filtrada.csv"  # Cambia la ruta al archivo CSV filtrado previamente

# Leer el archivo CSV con pandas
df = pd.read_csv(segudno_archivo_csv)

# Eliminar filas con "Standard Units" igual a "µg·mL⁻¹"
df_filtrado = df[df["Standard Units"] != "ug.mL-1"]

# Convertir "Standard Value" a tipo numérico y eliminar filas con valor igual a 0.0
df_filtrado["Standard Value"] = pd.to_numeric(df_filtrado["Standard Value"], errors='coerce')
df_filtrado = df_filtrado[df_filtrado["Standard Value"] != 0.0]

# Guardar el archivo CSV filtrado
dos_nuevo_csv = "archivo_filtrado.csv"
df_filtrado.to_csv(dos_nuevo_csv, index=False)

print(f"Archivo CSV filtrado generado: {dos_nuevo_csv}")
