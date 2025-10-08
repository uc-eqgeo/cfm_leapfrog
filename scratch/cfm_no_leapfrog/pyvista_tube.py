import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
# Create points for a line in the xy plane with a slight bend
n_points = 100
theta = np.linspace(0, np.pi / 4, n_points)  # Range from 0 to 45 degrees for a slight bend
x = np.linspace(0, 10, n_points)
y = np.sin(theta) * 4  # Add a slight bend
z = np.zeros_like(x)  # Keep z at zero (xy plane)

# Stack the points to create a line
line_points = np.column_stack((x, y, z))



# plot using pyplot
plt.plot(x, y)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()
# # Create a PyVista polyline
# line = pv.PolyData(line_points)
# line = line.cast_to_unstructured_grid()  # Cast to UnstructuredGrid

# # Visualize the line
# plotter = pv.Plotter()
# plotter.add_mesh(line, color='blue', line_width=5)
# plotter.show_grid()
# plotter.show()