import mysql.connector
import Credenciales
from mysql.connector import Error

try:
    connection = mysql.connector.connect(
        host='132.248.220.31',       # O la IP de tu servidor
        user=Credenciales.usuario,
        password=Credenciales.contrasena
    )

    if connection.is_connected():
        db_Info = connection.get_server_info()
        print(f"Conectado exitosamente a MySQL Server versión {db_Info}")
        
        # Aquí puedes ejecutar consultas, actualizaciones, creaciones o borrados
        cursor = connection.cursor()
        #cursor.execute("...")
        cursor.execute("SHOW DATABASES")

    for x in cursor:
          print(x)

except Error as e:
    print(f"Error al conectar a MySQL: {e}")

finally:
    if 'connection' in locals() and connection.is_connected():
        cursor.close()
        connection.close()
        print("La conexión a MySQL está cerrada.")
