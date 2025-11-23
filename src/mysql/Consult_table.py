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



print ("#################################")
print ("El contenido completo de la tabla")
cursor.execute("SELECT * FROM ExpDiff WHERE ExpDiff_id>7000")

myresult = cursor.fetchall()

for x in myresult:
  print(x)

print ("#################################")
print ("consultando algunas columnas de tabla")

cursor.execute("SELECT ExpDiff_id,gene_id,Log2FC FROM ExpDiff WHERE ExpDiff_id<10")

myresult = cursor.fetchall()

for x in myresult:
  print(x)

print ("#################################")
print ("consultando solo un registo de  tabla")
cursor.execute("SELECT * FROM ExpDiff WHERE gene_id='ID=FBgn0263392'")

myresult = cursor.fetchone()

print(myresult)

cursor.close()
conn.close()
