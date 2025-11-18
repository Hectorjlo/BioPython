import mysql.connector
import Credenciales

mydb = mysql.connector.connect(
  host=Credenciales.host,
  user=Credenciales.usuario,
  password=Credenciales.contrasena,
  database=Credenciales.database
)

# La sentencia SQL para crear la tabla ExpDiff
# Se añade "IF NOT EXISTS" para evitar errores si la tabla ya existe
CREATE_TABLE_QUERY = """
CREATE TABLE IF NOT EXISTS ExpDiff (
    ExpDiff_id INT AUTO_INCREMENT PRIMARY KEY,
    gene_id VARCHAR(50) UNIQUE,
    baseMean DECIMAL(40, 30),
    log2FC DECIMAL(40, 30),
    lfcSE DECIMAL(40, 30),
    stats DECIMAL(40, 30),
    pvalue DOUBLE(40,30),
    padj DOUBLE(40,30)
);
"""

mycursor = mydb.cursor()

mycursor.execute(CREATE_TABLE_QUERY)
print("La tabla 'ExpDiff' ha sido creada exitosamente o ya existía.")

if mycursor:
     mycursor.close()
     if mydb and mydb.is_connected():
            mydb.close()
            print("Conexión a MySQL cerrada.")
