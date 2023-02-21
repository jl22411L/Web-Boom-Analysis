# Web-Boom-Analysis of Single Cell Structures
This MATLAB code performs a structural analysis using the Web-Boom idealisation. 

By inputting the x and y positions of the booms, along with the loads that the structure experiences, you can find the shear flows in the Webs and direct stresses in the boom.
This Web-Boom Idealisation is usually performed in aircraft structures. Aircraft Structures for Engineering Students  by Megson, T.H.G. gives more information in the process.

This code (currently) only works on single cell structures. It also makes the assumption that the webs are straight lines between the booms. So you can't input curved surfaces for the webs.

The Booms are numberes clockwise. So if you have your booms layed out and pick a starting boom the next boom inputted will be the one connecting in the clockwise direction. Also, when you input the boom conditions the ith boom is connected to the i-1 and i+1 boom with web j-1 and j respectively. In short, make sure that the position of the booms makes sense.

Within the code, you can also make changes to the individual areas, A, and Youngs moduluses, E, of the booms, along with the thickness, t, and Shear Modulus, G, of the webs.

Note: when inputting the parameters for the different webs and booms, the row of the array corresponds to the same individual web or boom.

A = [Boom_1, Boom_2, ..., Boom_n]

E = [Boom_1, Boom_2, ..., Boom_n]


t = [Web_1, Web_2, ..., Web_n]

G = [Web_1, Web_2, ..., Web_n]


Further, their is another piece of code which can find the position of booms within a NACA 4 series aerofoil. By inputting the properties of the aerofoil such as max camber, position of max camber, and thickness to chord ratio it a function and then finds the y position that corresponds to the desired x positions of the booms.
