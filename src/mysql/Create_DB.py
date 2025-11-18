import mysql.connector
import Credenciales

# 1. Establecer la conexión a MySQL (sin especificar una base de datos aún)
# Asegúrate de reemplazar 'tu_usuario' y 'tu_contraseña' con tus credenciales
try:
    conn = mysql.connector.connect(
        host=Credenciales.host,
        user=Credenciales.usuario,
        password=Credenciales.contrasena
    )
    print("Conexión a MySQL exitosa.")

    # 2. Crear un objeto cursor
    cursor = conn.cursor()

    # 3. Ejecutar el comando SQL para crear la base de datos
    nombre_db = "pasilla"
    cursor.execute(f"CREATE DATABASE {nombre_db}")
    print(f"Base de datos '{nombre_db}' creada exitosamente.")

except mysql.connector.Error as error:
    print(f"Error al crear la base de datos: {error}")

finally:
    # 4. Cerrar el cursor y la conexión
    if cursor:
        cursor.close()
    if conn and conn.is_connected():
        conn.close()
        print("La conexión a MySQL se ha cerrado.")
