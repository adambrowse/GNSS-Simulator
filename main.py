import sys
import numpy as np
import math
import random
from PyQt5.QtWidgets import QSlider, QFrame, QScrollArea, QMessageBox, QDoubleSpinBox, QLineEdit, QLabel, QDialog, QPushButton, QCheckBox, QApplication, QMainWindow, QOpenGLWidget, QHBoxLayout, QVBoxLayout, QWidget
from PyQt5.QtCore import QTimer, Qt
from OpenGL.GL import *
from OpenGL.GLU import *


class Camera:
    def __init__(self):
        self.distance = 70000
        self.azimuth = 0
        self.elevation = -60
        self.skew = 0
        self.aspect_ratio = 4/3

    def apply(self, aspect_ratio=None):
        if aspect_ratio is not None:
            self.aspect_ratio = aspect_ratio

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, self.aspect_ratio, 1, 1000000)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef(0, 0, -self.distance)
        glRotatef(self.elevation, 1, 0, 0)
        glRotatef(self.azimuth, 0, 1, 0)
        glRotatef(self.skew, 0, 0, 1)

    def move_forward(self, factor):
        self.distance -= self.distance*factor

    def move_backward(self, factor):
        self.distance += self.distance*factor

    def rotate_x_plus(self, angle):
        self.elevation -= angle

    def rotate_x_minus(self, angle):
        self.elevation += angle

    def rotate_z_plus(self, angle):
        self.skew += angle

    def rotate_z_minus(self, angle):
        self.skew -= angle


class gnssSimulator(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.camera = Camera()

        # Define physical parameters
        self.earth_radius = 6371
        self.gravitational_parameter = 398600

        # Create timer for animation
        self.speed_up = 30
        self.timer = QTimer()
        self.timer.timeout.connect(self.animate_satellite)
        self.timer.start(30)
        self.timer.stop()

        # Define lat and long of receiver
        self.latitude = 30
        self.longitude = -90

        # Initialise the list of satellites
        self.satellites = [
            {"name": "Sat1", "true_anomaly": math.radians(0), "eccentricity": 0, "semi_major_axis": 10000, "inclination": math.radians(55), "raan": math.radians(0), "colour": (1, 0, 0)}, 
            {"name": "Sat2", "true_anomaly": math.radians(60), "eccentricity": 0, "semi_major_axis": 10000, "inclination": math.radians(55), "raan": math.radians(60), "colour": (0, 1, 0)},
            {"name": "Sat3", "true_anomaly": math.radians(120), "eccentricity": 0, "semi_major_axis": 10000, "inclination": math.radians(55), "raan": math.radians(120), "colour": (1, 0.7, 0)}, 
            {"name": "Sat4", "true_anomaly": math.radians(180), "eccentricity": 0, "semi_major_axis": 10000, "inclination": math.radians(55), "raan": math.radians(180), "colour": (1, 1, 0)},
            {"name": "Sat5", "true_anomaly": math.radians(240), "eccentricity": 0, "semi_major_axis": 10000, "inclination": math.radians(55), "raan": math.radians(240), "colour": (1, 0, 1)}, 
            {"name": "Sat6", "true_anomaly": math.radians(300), "eccentricity": 0, "semi_major_axis": 10000, "inclination": math.radians(55), "raan": math.radians(300), "colour": (0, 1, 1)}, 
        ]

        # State acceptable satellite colours
        self.suitable_colours = [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            (1, 1, 0),
            (1, 0, 1),
            (0, 1, 1)
        ]

        # Set visibility of objects
        self.receiver_on = False
        self.satellites_on = False
        self.orbits_on = False
        self.spheres_on = False
        self.lines_of_sight_on = False
        self.circles_on = False

        # Add toggle buttons
        self.toggle_layout = QVBoxLayout()
        self.add_toggle_button('Orbits', self.toggle_orbits)
        self.add_toggle_button('Satellites', self.toggle_satellites)
        self.add_toggle_button('Receiver', self.toggle_receiver)
        self.add_toggle_button('Spheres', self.toggle_spheres)
        self.add_toggle_button('Lines of sight', self.toggle_lines_of_sight)
        self.add_toggle_button('Circles', self.toggle_circles)

        # Add a slider to control the animation speed
        self.speed_slider = QSlider(Qt.Horizontal)
        self.speed_slider.setRange(1, 300)
        self.speed_slider.setValue(self.speed_up)
        self.speed_slider.valueChanged.connect(self.update_speed_up)
        self.toggle_layout.addWidget(self.speed_slider)

        # Create a QWidget to contain the toggles
        self.toggle_widget = QWidget()
        self.toggle_widget.setLayout(self.toggle_layout)

        # Add the satellite editing buttons
        self.button_layout = QVBoxLayout()
        self.satellite_form_buttons()

        # Add the widgets
        main_layout = QVBoxLayout()
        main_layout.addWidget(self.toggle_widget)
        self.setLayout(main_layout)

    def animate_satellite(self):
        for orbit_params in self.satellites:
            # Update satellite angle for animation
            v = orbit_params["true_anomaly"]
            e = orbit_params["eccentricity"]
            a = orbit_params["semi_major_axis"]
            mu = self.gravitational_parameter

            r = a * (1 - e ** 2) / (1 + e * math.cos(v))
            omega = math.sqrt(mu / r ** 3)

            v += omega * self.speed_up

            if v >= 2 * np.pi:
                v = 0

            orbit_params["true_anomaly"] = v

        self.update()

    def draw_lines_of_latitude(self):
        radius = self.earth_radius + 80

        glColor3f(0, 0, 1)

        # Draw the latitude lines every 10 degrees
        for i in range(-90, 91, 10):
            angle = math.radians(i)
            
            for j in range(0, 360, 10):
                glBegin(GL_LINE_STRIP)
                x1 = radius * math.cos(math.radians(j)) * math.cos(angle)
                y1 = radius * math.sin(math.radians(j)) * math.cos(angle)
                z1 = radius * math.sin(angle)
                x2 = radius * math.cos(math.radians(j + 10)) * math.cos(angle)
                y2 = radius * math.sin(math.radians(j + 10)) * math.cos(angle)
                z2 = radius * math.sin(angle)
                glVertex3f(x1, y1, z1)
                glVertex3f(x2, y2, z2)
                glEnd()

    def draw_lines_of_longitude(self):
        radius = self.earth_radius + 80

        glColor3f(0, 0, 1)

        # Draw the longitude lines every 10 degrees
        for i in range(-180, 181, 10): 
            angle = math.radians(i)
            
            for j in range(-90, 91, 10):
                glBegin(GL_LINE_STRIP)
                x1 = radius * math.cos(math.radians(j)) * math.cos(angle)
                y1 = radius * math.sin(math.radians(j)) * math.cos(angle)
                z1 = radius * math.sin(angle)
                x2 = radius * math.cos(math.radians(j)) * math.cos(angle + math.radians(10))
                y2 = radius * math.sin(math.radians(j)) * math.cos(angle + math.radians(10))
                z2 = radius * math.sin(angle + math.radians(10))
                glVertex3f(x1, y1, z1)
                glVertex3f(x2, y2, z2)
                glEnd()

    def initializeGL(self):
        # Initialise the lighting
        glClearColor(0, 0, 0, 1)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE)

        # Enable lighting
        glEnable(GL_LIGHTING)

        # Set the light position (above-right-front)
        light_position = [1, 1, 1, 0]
        glLightfv(GL_LIGHT0, GL_POSITION, light_position)
        light_diffuse = [1, 1, 1, 1]
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse)

        glEnable(GL_LIGHT0)

        # Set the lighting for the Earth
        earth_ambient = [0.2, 0.2, 0.5, 1]
        earth_diffuse = [0.2, 0.2, 0.6, 1]
        earth_specular = [0.1, 0.1, 0.1, 1]
        earth_shininess = 10
        glMaterialfv(GL_FRONT, GL_AMBIENT, earth_ambient)
        glMaterialfv(GL_FRONT, GL_DIFFUSE, earth_diffuse)
        glMaterialfv(GL_FRONT, GL_SPECULAR, earth_specular)
        glMaterialfv(GL_FRONT, GL_SHININESS, earth_shininess)

    def paintGL(self):
        try:
            # Style components
            for i in range(self.toggle_layout.count()):
                toggle_button = self.toggle_layout.itemAt(i).widget()
                toggle_button.setStyleSheet("""
                    QCheckBox::indicator {
                        width: 40px;
                        height: 20px;
                        border-radius: 10px;
                        background-color: #ddd;
                    }
                    QCheckBox::indicator:checked {
                        background-color: #00BFFF;
                    }
                    QCheckBox {
                        color: white;
                    }
                """)
            
            super().paintGL()
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            self.camera.apply()

            # Draw earth
            self.draw_earth()

            # Draw receiver
            if self.receiver_on:
                x_rec, y_rec, z_rec = self.draw_receiver()

            # Draw the satellite orbit
            if self.orbits_on:
                for orbit_params in self.satellites:
                    self.draw_satellite_orbit(orbit_params)

            # Empty array for satellite coordinates
            sat_coords = []

            # Draw the satellites, receiver, lines of sight and circles
            for orbit_params in self.satellites:
                if self.satellites_on:
                    x_sat, y_sat, z_sat = self.draw_satellite(orbit_params)
                    sat_coords.append([x_sat, y_sat, z_sat])

                    # Check if the receiver exists and has an uninterrupted view of satellite
                    if self.receiver_on and self.has_uninterrupted_view((x_sat, y_sat, z_sat), (x_rec, y_rec, z_rec)):
                        if self.lines_of_sight_on:
                            self.draw_connection_line((x_sat, y_sat, z_sat), (x_rec, y_rec, z_rec), orbit_params)
                        if self.circles_on:
                            self.draw_circle_on_earth(orbit_params, (x_sat, y_sat, z_sat), (x_rec, y_rec, z_rec))
                            
            # Draw the spheres
            for i, orbit_params in enumerate(self.satellites):
                if self.satellites_on:
                    x_sat, y_sat, z_sat = sat_coords[i]

                    # Check if the receiver exists and has an uninterrupted view of satellite
                    if self.receiver_on and self.has_uninterrupted_view((x_sat, y_sat, z_sat), (x_rec, y_rec, z_rec)):
                        if self.spheres_on:
                            self.draw_satellite_range(orbit_params, (x_sat, y_sat, z_sat), (x_rec, y_rec, z_rec))

            # Draw the toggle widget
            self.toggle_widget.resize(200, 200)
            self.toggle_widget.move(10, self.height() - 210)
            self.toggle_widget.show()
        
        except Exception as e:
            QMessageBox.warning(self, "Error", str(e))

    def draw_earth(self):
        # Translate the Earth to the center
        glTranslatef(0, 0, 0)

        glEnable(GL_LIGHTING)

        # Draw the Earth
        glColor3f(0, 0, 1)
        quadric = gluNewQuadric()
        gluSphere(quadric, self.earth_radius, 100, 100)

        glDisable(GL_LIGHTING)

        # Draw the lines of latitude and longitude
        self.draw_lines_of_latitude()
        self.draw_lines_of_longitude()

    def draw_receiver(self):
        # Calculate the position of the receiver
        latitude = self.latitude
        longitude = self.longitude
        radius = self.earth_radius + 50
        x = radius * math.cos(math.radians(longitude)) * math.cos(math.radians(latitude))
        y = radius * math.sin(math.radians(longitude)) * math.cos(math.radians(latitude))
        z = radius * math.sin(math.radians(latitude))

        glTranslatef(x, y, z)

        glColor3f(1, 0, 0)
        quadric = gluNewQuadric()
        gluSphere(quadric, 200, 100, 100)

        # Translate back to the origin
        glTranslatef(-x, -y, -z)

        return(x, y, z)

    def draw_satellite(self, orbit_params):
        # Calculate satellite's position based on its current angle
        v = orbit_params["true_anomaly"]
        a = orbit_params["semi_major_axis"]
        e = orbit_params["eccentricity"]
        i = orbit_params["inclination"]
        raan = orbit_params["raan"]
        colour = orbit_params["colour"]
        r = a * (1 - e ** 2) / (1 + e * math.cos(v))
        x = r * (math.cos(raan) * math.cos(v) - math.sin(raan) * math.sin(v) * math.cos(i))
        y = r * (math.sin(raan) * math.cos(v) + math.cos(raan) * math.sin(v) * math.cos(i))
        z = r * math.sin(v) * math.sin(i)

        # Draw the satellite
        glTranslatef(x, y, z)
        R = 200
        glColor3f(*colour)
        quadric = gluNewQuadric()
        gluSphere(quadric, R, 100, 100)

        glTranslatef(-x, -y, -z)

        return(x, y, z)

    def draw_circle_on_earth(self, orbit_params, satellite_position, receiver_position):
        # Convert points to numpy arrays
        origin = np.array((0, 0, 0))
        receiver = np.array(receiver_position)
        satellite = np.array(satellite_position)

        # Calculate direction vector from origin to satellite
        normal_vector = satellite - origin

        # Normalise the normal vector
        normal_vector = normal_vector / np.linalg.norm(normal_vector)

        # Calculate the vector from the origin to the receiver
        vector_to_receiver = receiver - origin

        # Calculate the projection of vector_to_receiver on the direction vector
        scalar_projection = np.dot(vector_to_receiver, normal_vector) / np.dot(normal_vector, normal_vector)
        closest_point = origin + scalar_projection * normal_vector

        radius = np.linalg.norm(np.array(closest_point) - np.array(receiver_position))

        A, B, C = normal_vector
        D = -np.dot(normal_vector, closest_point)

        # Set the number of points to create the circle
        num_points = 100
        
        # Create the vertex array
        vertex_array = np.empty((num_points + 1, 3), dtype=np.float32)
        
        if B == 0 and C == 0:
            basis_vector1 = np.array([1, 0, 0])
        else:
            basis_vector1 = np.array([-B, A, 0])

        # Find basis vectors
        basis_vector1 = basis_vector1 / np.linalg.norm(basis_vector1)
        basis_vector2 = np.cross(normal_vector, basis_vector1)
        basis_vector2 = basis_vector2 / np.linalg.norm(basis_vector2)

        # Create the array of points for the circle
        for i in range(num_points + 1):
            angle = i * 2 * np.pi / num_points
            x = closest_point[0] + radius * (basis_vector1[0] * np.cos(angle) + basis_vector2[0] * np.sin(angle))
            y = closest_point[1] + radius * (basis_vector1[1] * np.cos(angle) + basis_vector2[1] * np.sin(angle))
            z = closest_point[2] + radius * (basis_vector1[2] * np.cos(angle) + basis_vector2[2] * np.sin(angle))
            vertex_array[i] = [x, y, z]

        colour = orbit_params["colour"]
        glColor3f(*colour)
        
        # Draw the circle
        glBegin(GL_LINE_LOOP)
        
        for i in range(num_points + 1):
            glVertex3fv(vertex_array[i])

        glEnd()

    def draw_satellite_range(self, orbit_params, satellite_position, receiver_position):
        # Get cartesian positions
        receiver = np.array(receiver_position)
        satellite = np.array(satellite_position)
        radius = np.linalg.norm(satellite - receiver)

        # Calculate satellite's position based on its current angle
        v = orbit_params["true_anomaly"]
        a = orbit_params["semi_major_axis"]
        e = orbit_params["eccentricity"]
        i = orbit_params["inclination"]
        raan = orbit_params["raan"]
        colour = orbit_params["colour"]
        r = a * (1 - e ** 2) / (1 + e * math.cos(v))
        x = r * (math.cos(raan) * math.cos(v) - math.sin(raan) * math.sin(v) * math.cos(i))
        y = r * (math.sin(raan) * math.cos(v) + math.cos(raan) * math.sin(v) * math.cos(i))
        z = r * math.sin(v) * math.sin(i)

        # Set satellite's position
        glTranslatef(x, y, z)
        R = radius
        glColor4f(*colour, 0.25)

        # Draw a semi-opaque sphere
        quadric = gluNewQuadric()
        gluSphere(quadric, R, 100, 100)

        glTranslatef(-x, -y, -z)

    def draw_satellite_orbit(self, orbit_params):
        # Draw satellite orbit
        glColor3f(1, 1, 1)

        glBegin(GL_LINE_LOOP)
        for ang in range(360):
            v = math.radians(ang)
            a = orbit_params["semi_major_axis"]
            e = orbit_params["eccentricity"]
            i = orbit_params["inclination"]
            raan = orbit_params["raan"]
            r = a * (1 - e ** 2) / (1 + e * math.cos(v))
            x = r * (math.cos(raan) * math.cos(v) - math.sin(raan) * math.sin(v) * math.cos(i))
            y = r * (math.sin(raan) * math.cos(v) + math.cos(raan) * math.sin(v) * math.cos(i))
            z = r * math.sin(v) * math.sin(i)

            glVertex3f(x, y, z)

        glEnd()

    def draw_connection_line(self, satellite_position, receiver_position, orbit_params):
        # Draw the line
        colour = orbit_params["colour"]
        glColor3f(*colour)
        glBegin(GL_LINES)
        glVertex3f(satellite_position[0], satellite_position[1], satellite_position[2])
        glVertex3f(receiver_position[0], receiver_position[1], receiver_position[2])
        glEnd()

    def resizeGL(self, w, h):
        # Set the viewport
        glViewport(0, 0, w, h)

        aspect_ratio = w / h

        # Set the projection matrix
        self.camera.apply(aspect_ratio)

    def has_uninterrupted_view(self, satellite_position, receiver_position):
        # Define numpy arrays for positions
        satellite_position = np.array(satellite_position) 
        receiver_position = np.array(receiver_position)

        # Set the direction vector between the points
        direction_vector = receiver_position - satellite_position

        # Move from the satellite to the receiver in small steps, if distance to origin is less than the earth's radius then no connection 
        steps = 100
        for i in range(steps):
            point = satellite_position + (i/steps)*direction_vector
            if np.linalg.norm(point) < self.earth_radius:
                return False
            else:
                pass

        return True

    def wheelEvent(self, event):
        # Define zoom in/out using mouse wheel
        if event.angleDelta().y() > 0:
            self.camera.move_forward(0.1)
            self.update()

        elif event.angleDelta().y() < 0:
            self.camera.move_backward(0.1)
            self.update()

    def keyPressEvent(self, event):
        # Move camera using keys
        if event.key() == Qt.Key_D:
            self.camera.rotate_z_plus(2)
            self.update()
        elif event.key() == Qt.Key_A:
            self.camera.rotate_z_minus(2)
            self.update()
        elif event.key() == Qt.Key_W:
            self.camera.rotate_x_plus(5)
            self.update()
        elif event.key() == Qt.Key_S:
            self.camera.rotate_x_minus(5)
            self.update()
        elif event.key() == Qt.Key_Space:
            if self.timer.isActive():
                self.timer.stop()
            else:
                self.timer.start()

    def add_toggle_button(self, text, method):
        # Create a new toggle button
        toggle_button = QCheckBox(text)
        
        # Ignore keyboard presses eg space bar
        toggle_button.setFocusPolicy(Qt.NoFocus)

        # Connect the button to a function
        toggle_button.stateChanged.connect(method)

        # Add the button to the layout
        self.toggle_layout.addWidget(toggle_button)

        # Set the initial state of the button based on its text
        if text in ["Satellites", "Orbits"]:
            toggle_button.setChecked(True)
        else:
            toggle_button.setChecked(False)

    def toggle_orbits(self):
        self.orbits_on = not self.orbits_on
        self.update()

    def toggle_satellites(self):
        self.satellites_on = not self.satellites_on
        self.update()

    def toggle_receiver(self):
        self.receiver_on = not self.receiver_on
        self.update()

    def toggle_spheres(self):
        self.spheres_on = not self.spheres_on
        self.update()

    def toggle_lines_of_sight(self):
        self.lines_of_sight_on = not self.lines_of_sight_on
        self.update()

    def toggle_circles(self):
        self.circles_on = not self.circles_on
        self.update()

    def satellite_form_buttons(self):
        self.open_form_button = QPushButton("Add Satellite", self)
        self.open_form_button.setGeometry(10, 10, 100, 30)
        self.open_form_button.clicked.connect(self.new_satellite_form)

        # Ignore keyboard presses eg space bar
        self.open_form_button.setFocusPolicy(Qt.NoFocus)

        self.view_satellites_button = QPushButton("View Satellites", self)
        self.view_satellites_button.setGeometry(10, 50, 100, 30)
        self.view_satellites_button.clicked.connect(self.view_satellites_form)

        # Ignore keyboard presses eg space bar
        self.view_satellites_button.setFocusPolicy(Qt.NoFocus)

        self.edit_receiver_button = QPushButton("Edit Receiver", self)
        self.edit_receiver_button.setGeometry(10, 90, 100, 30)
        self.edit_receiver_button.clicked.connect(self.edit_receiver_form)

        # Ignore keyboard presses eg space bar
        self.edit_receiver_button.setFocusPolicy(Qt.NoFocus)

    def new_satellite_form(self):
        # Create a new window
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Satellite")
        dialog.setGeometry(200, 200, 300, 200)
        layout = QVBoxLayout()

        # Create input fields for each satellite parameter
        name_label = QLabel("Name:")
        name_input = QLineEdit()
        name_input.setText("New Satellite")
        layout.addWidget(name_label)
        layout.addWidget(name_input)

        true_anomaly_label = QLabel("True Anomaly (degrees):")
        true_anomaly_input = QDoubleSpinBox()
        true_anomaly_input.setRange(0, 360)
        true_anomaly_input.setSingleStep(1)
        layout.addWidget(true_anomaly_label)
        layout.addWidget(true_anomaly_input)

        eccentricity_label = QLabel("Eccentricity:")
        eccentricity_input = QDoubleSpinBox()
        eccentricity_input.setRange(0, 1)
        eccentricity_input.setSingleStep(0.01)
        layout.addWidget(eccentricity_label)
        layout.addWidget(eccentricity_input)

        semi_major_axis_label = QLabel("Semi-Major Axis (km):")
        semi_major_axis_input = QDoubleSpinBox()
        semi_major_axis_input.setRange(6700, 1000000) 
        semi_major_axis_input.setSingleStep(100)
        layout.addWidget(semi_major_axis_label)
        layout.addWidget(semi_major_axis_input)

        inclination_label = QLabel("Inclination (degrees):")
        inclination_input = QDoubleSpinBox()
        inclination_input.setRange(0, 360)
        inclination_input.setSingleStep(1)
        layout.addWidget(inclination_label)
        layout.addWidget(inclination_input)

        raan_label = QLabel("RAAN (degrees):")
        raan_input = QDoubleSpinBox()
        raan_input.setRange(0, 360)
        raan_input.setSingleStep(1)
        layout.addWidget(raan_label)
        layout.addWidget(raan_input)

        # Create an "OK" button to confirm the input values
        ok_button = QPushButton("OK")
        def handle_ok_button_click():
            for satellite in self.satellites:
                if satellite["name"] == name_input.text():
                    QMessageBox.warning(None, "Choose another name", "This satellite name already exists")
                    return
                
            colour = self.suitable_colours[random.randint(0, 5)]
            self.add_satellite_from_form(
                name_input.text(),
                true_anomaly_input.value(),
                eccentricity_input.value(),
                semi_major_axis_input.value(),
                inclination_input.value(),
                raan_input.value(),
                colour
            )
            dialog.accept()
        ok_button.clicked.connect(handle_ok_button_click)
        layout.addWidget(ok_button)

        dialog.setLayout(layout)

        dialog.exec_()

    def add_satellite_from_form(self, name, true_anomaly, eccentricity, semi_major_axis, inclination, raan, colour):
        # Convert true anomaly and inclination from degrees to radians
        true_anomaly = math.radians(true_anomaly)
        inclination = math.radians(inclination)
        raan = math.radians(raan)

        if eccentricity == 1:
            eccentricity = 0.99

        # Add the satellite to the list
        self.satellites.append({"name": name, 
                                "true_anomaly": true_anomaly, 
                                "eccentricity": eccentricity, 
                                "semi_major_axis": semi_major_axis, 
                                "inclination": inclination, 
                                "raan": raan, 
                                "colour": colour})

        self.update()

    def view_satellites_form(self):
        # Pause the animation if it's going
        if self.timer.isActive():
            self.timer.stop()

        # Create a new window
        dialog = QDialog(self)
        dialog.setWindowTitle("View Satellites")
        dialog.setGeometry(200, 200, 400, 400)
        layout = QVBoxLayout()

        label = QLabel("Name, Distance to receiver (km), colour")
        label.setAlignment(Qt.AlignLeft)
        layout.addWidget(label)

        # Area to hold the list of satellites
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout()
        scroll_widget.setLayout(scroll_layout)
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area)

        # Widget to hold the satellites
        satellites_widget = QWidget()
        satellites_layout = QVBoxLayout()
        satellites_widget.setLayout(satellites_layout)

        # Function for when a satellite is deleted
        def delete_satellite(index, label, button, box, distance):
            index = satellite_label_list.index(label)
            satellite_label_list.pop(index)
            
            del self.satellites[index]

            label.deleteLater()
            button.deleteLater()
            box.deleteLater()
            distance.deleteLater()

            self.update()

        satellite_label_list = []

        receiver_position = self.find_receiver_position()

        for index, satellite in enumerate(self.satellites):
            satellite_name = satellite["name"]
            satellite_label = QLabel(satellite_name)
            delete_button = QPushButton("Delete")

            colour_box = QFrame()
            colour_box.setFrameShape(QFrame.Box)
            colour_box.setFixedSize(20, 20)
            colour = satellite["colour"]
            R, G, B = int(colour[0]), int(colour[1]), int(colour[2])
            max_val = max(R, G, B)
            if max_val == 0:
                hex_colour = '#000000'
            else:
                R_norm = R / max_val * 255
                G_norm = G / max_val * 255
                B_norm = B / max_val * 255
                hex_colour = '#{:02x}{:02x}{:02x}'.format(int(R_norm), int(G_norm), int(B_norm))
            colour_box.setStyleSheet("background-color: {}".format(hex_colour))

            satellite_label_list.append(satellite_label)

            satellite_position = self.find_satellite_position(satellite)
            
            # Only show distance is the satellite can see the receiver
            if self.receiver_on and self.has_uninterrupted_view(satellite_position, receiver_position):
                distance_label = QLabel(self.find_distance(satellite_position, receiver_position))
            else:
                distance_label = QLabel("No connection")

            delete_button.clicked.connect(lambda _, 
                                          i=index, 
                                          lab=satellite_label, 
                                          btn=delete_button, 
                                          box=colour_box, 
                                          dis=distance_label: 
                                          delete_satellite(i, lab, btn, box, dis))
            
            satellite_layout = QHBoxLayout()
            satellite_layout.addWidget(satellite_label)
            satellite_layout.addWidget(distance_label)
            satellite_layout.addWidget(colour_box)
            satellite_layout.addWidget(delete_button)
            satellites_layout.addLayout(satellite_layout)

        satellites_layout.addStretch()
        scroll_layout.addWidget(satellites_widget)

        # Create an "OK" button to close the dialog
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(dialog.accept)
        layout.addWidget(ok_button)

        dialog.setLayout(layout)

        dialog.exec_()

    def edit_receiver_form(self):
        # Create a new window
        dialog = QDialog(self)
        dialog.setWindowTitle("Edit Receiver")
        dialog.setGeometry(200, 200, 300, 200)
        layout = QVBoxLayout()

        # Create input fields for latitude and longitude
        latitude_label = QLabel("Latitude:")
        latitude_input = QDoubleSpinBox()
        latitude_input.setRange(-90, 90)
        latitude_input.setSingleStep(1)
        latitude_input.setValue(self.latitude)
        layout.addWidget(latitude_label)
        layout.addWidget(latitude_input)

        longitude_label = QLabel("Longitude:")
        longitude_input = QDoubleSpinBox()
        longitude_input.setRange(-180, 180)
        longitude_input.setSingleStep(1)
        longitude_input.setValue(self.longitude)
        layout.addWidget(longitude_label)
        layout.addWidget(longitude_input)

        # Create an "OK" button to confirm the input values
        ok_button = QPushButton("OK")
        def handle_ok_button_click():
            self.latitude = latitude_input.value()
            self.longitude = longitude_input.value()
            dialog.accept()
        ok_button.clicked.connect(handle_ok_button_click)
        layout.addWidget(ok_button)

        dialog.setLayout(layout)

        dialog.exec_()

    def update_speed_up(self, value):
        self.speed_up = value

    def find_receiver_position(self):
        # Calculate the position of the receiver
        latitude = self.latitude
        longitude = self.longitude
        radius = self.earth_radius + 10
        x_rec = radius * math.cos(math.radians(longitude)) * math.cos(math.radians(latitude))
        y_rec = radius * math.sin(math.radians(longitude)) * math.cos(math.radians(latitude))
        z_rec = radius * math.sin(math.radians(latitude))

        return np.array([x_rec, y_rec, z_rec])

    def find_satellite_position(self, satellite):
        # Calculate the position of the satellite
        v = satellite["true_anomaly"]
        a = satellite["semi_major_axis"]
        e = satellite["eccentricity"]
        i = satellite["inclination"]
        raan = satellite["raan"]
        r = a * (1 - e ** 2) / (1 + e * math.cos(v))
        x_sat = r * (math.cos(raan) * math.cos(v) - math.sin(raan) * math.sin(v) * math.cos(i))
        y_sat = r * (math.sin(raan) * math.cos(v) + math.cos(raan) * math.sin(v) * math.cos(i))
        z_sat = r * math.sin(v) * math.sin(i)

        return np.array([x_sat, y_sat, z_sat])

    def find_distance(self, satellite_position, receiver_position):
        return str(round(np.linalg.norm(satellite_position - receiver_position), 2))
 
 
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle('GNSS Visualiser')
        
        self.setGeometry(100, 100, 1600, 800)
        
        # Create the app
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout(self.central_widget)
        self.gnss_simulator = gnssSimulator(self.central_widget)
        self.layout.addWidget(self.gnss_simulator)

        # Key press event connection
        self.keyPressEvent = self.gnss_simulator.keyPressEvent
        self.wheelEvent = self.gnss_simulator.wheelEvent

# Call the app
app = QApplication(sys.argv)

# Create the main window
main_window = MainWindow()
main_window.show()

# Let's go!
sys.exit(app.exec_())
