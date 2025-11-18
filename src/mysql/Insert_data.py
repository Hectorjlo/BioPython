import mysql.connector
import csv
import Credenciales

conn = mysql.connector.connect(
  host=Credenciales.host,
  user=Credenciales.usuario,
  password=Credenciales.contrasena,
  database=Credenciales.database
)

cursor = conn.cursor()

# Ruta del archivo
archivo = '/export/space3/users/vjimenez/Data/ExpDifRes_v2.txt'
contador=0;

# Leer y procesar el archivo
with open(archivo, 'r') as file:
    file.readline() # Esto lee y descarta la primera línea
    for linea in file:
        valores = linea.strip().split(',')
        query = "INSERT INTO ExpDiff (gene_id,baseMean, log2FC, lfcSE, stats, pvalue, padj) VALUES (%s, %s,%s, %s, %s, %s, %s)"
        cursor.execute(query, valores) 
        contador= contador+cursor.rowcount

# Confirmar los cambios.
conn.commit()
# Imprimir el número de filas afectadas
print(contador, "registros insertados.")

cursor.close()
conn.close()
print ( "Done")

