import yaml
import numpy as np


def yaml_loader(filepath):
    """Reads in yaml files"""
    with open(filepath,"r") as file_descriptor:
        data = yaml.load(file_descriptor)
    return data

def load_frequencies(data):
    """Extracts a particular property from the disctionary and returns as a list"""
    new_list = np.zeros(shape=(len(data['phonon']),len(data['phonon'][0]['band'])))
    for i in range(len(data['phonon'])):
        for j in range(len(data['phonon'][i]['band'])):
            new_list[i,j] = float(data['phonon'][i]['band'][j]['frequency'])
    return new_list

def load_weights(data):
    """Extracts a particular property from the disctionary and returns as a list"""
    new_list = np.zeros(shape=(len(data['phonon'])))
    for i in range(len(data['phonon'])):
        new_list[i] = data['phonon'][i]['weight']
    return new_list

def load_positions(data):
    """Extracts a particular property from the disctionary and returns as a list"""
    positions = []
    for i in range(len(data['atoms'])):
	positions.append([data['atoms'][i]['position'],data['atoms'][i]['symbol'],data['atoms'][i]['mass']])
    return positions

def load_lattice(data):
    lattice = np.zeros(shape=(3,3))
    for i, element in enumerate(data['lattice']):
	lattice[i] = element
    return lattice

def load_q_points(data):
    """Extracts a particular property from the disctionary and returns as a list"""
    new_list = []
    for i in range(len(data['phonon'])):
        new_list.append(data['phonon'][i]['q-position'])
    return new_list

def load_eigenvectors(data):
    """Extracts a particular property from the disctionary and returns as a list"""
    w, h = len(data['phonon']), len(data['phonon'][0]['band'])
    w1 = len(data['phonon'][0]['band'][0]['eigenvector'])
    real_list = np.zeros(shape=(w,h,w1,3))
    im_list = np.zeros(shape=(w,h,w1,3))
    for i in range(w):
        for j in range(h):
             for k in range(w1):
                 for l in range(3):
                    real_list[i,j,k,l] = float(data['phonon'][i]['band'][j]['eigenvector'][k][l][0])
                    im_list[i,j,k,l] = data['phonon'][i]['band'][j]['eigenvector'][k][l][1]
    return real_list,im_list

data = yaml_loader('mesh.yaml')
# Load up the frequencies, they are stored a list of lists, sorted by q-point
frequencies = load_frequencies(data)
# Loads up the eigenvotors, they are sorted by [q-pt,band,atom,cartesian_direction]
real_vec, im_vec = load_eigenvectors(data)
# Loads up the q-points
q_points = load_q_points(data)
# Loads up the q-point weights
weights = load_weights(data)

positions = load_positions(data) 
lattice = load_lattice(data)

output = open('climax_input.phonon', 'w')
output.write(" BEGIN header \n")
output.write(" Number of ions         %i\n" % len(positions))
output.write(" Number of branches     %i\n" % len(data['phonon'][0]['band']))
output.write(" Number of wavevectors  %i\n" % len(data['phonon']))
output.write(" Frequencies in         cm-1\n")
output.write(" IR intensities in      (D/A)**2/amu\n")
output.write(" Raman intensities in   A**4\n")
output.write(" Unit cell vectors (A)\n")
for vector in lattice:
    output.write("   %8.5f   %8.5f   %8.5f \n" % (vector[0],vector[1],vector[2]))
output.write(" Fractional Co-ordinates\n")
for i, ion in enumerate(positions):
    output.write("  %i   %8.5f   %8.5f   %8.5f   %s   %8.5f \n" % (i+1,positions[i][0][0],
                positions[i][0][1],positions[i][0][2], positions[i][1], positions[i][2]))
output.write(" END header \n")

for i in range(len(data['phonon'])):
    output.write("    q-pt=   %i   %8.5f   %8.5f   %8.5f   %8.5f \n  " % ( i+1, q_points[i][0], q_points[i][1],
                                                                          q_points[i][2], weights[i]))
    for j in range(len(data['phonon'][0]['band'])):
        output.write("       %i   %8.5f \n" % (j, frequencies[i,j] * 33.35641))
    output.write("                  Phonon Eigenvectors   \n")
    output.write("Mode Ion                X                                   Y                                   Z\n")
    for j in range(len(data['phonon'][0]['band'])):
        for k in range(len(positions)):
            output.write("  %i  %i  %10.6f %10.6f   %10.6f %10.6f   %10.6f %10.6f\n" % (j+1, k+1, real_vec[i,j,k,0],
                        im_vec[i,j,k,0],real_vec[i,j,k,1],im_vec[i,j,k,1],real_vec[i,j,k,2],im_vec[i,j,k,2]))
output.close()
