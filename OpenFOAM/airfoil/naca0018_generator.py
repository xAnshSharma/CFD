import numpy as np
import math

def generate_bent_airfoil_coordinates(
    n_points=100, 
    chord_length=1.0, 
    x_offset=0.0, 
    y_offset=0.0, 
    angle_of_attack=0.0,
    bend_factor=0.0,
    thickness_front=0.18, # New parameter for thickness at the leading edge
    thickness_end=0.05    # New parameter for thickness at the trailing edge
):
    """
    Generate NACA0018 airfoil coordinates and apply a bending function.
    
    The bending is created by applying a parabolic camber line to the original airfoil shape.
    
    Parameters:
    - n_points: Number of points to define the airfoil.
    - chord_length: The length of the airfoil's chord.
    - x_offset, y_offset: Translation of the airfoil's origin.
    - angle_of_attack: Rotation of the airfoil in degrees.
    - bend_factor: A float that determines the amount of parabolic bending.
                   A value of 0.0 results in a straight airfoil.
    - thickness_front: The maximum thickness at the leading edge.
    - thickness_end: The thickness at the trailing edge.
    """
    # Use cosine spacing for better resolution near the leading edge
    beta = np.linspace(0, np.pi, n_points)
    x = 0.5 * chord_length * (1 - np.cos(beta))
    
    # Calculate a linear taper for the thickness
    thickness_taper = thickness_front * (1 - x / chord_length) + thickness_end * (x / chord_length)

    # NACA0018 thickness distribution with a variable thickness
    y_t = 5 * thickness_taper * chord_length * (
        0.2969 * np.sqrt(x / chord_length) -
        0.1260 * (x / chord_length) -
        0.3516 * (x / chord_length)**2 +
        0.2843 * (x / chord_length)**3 -
        0.1015 * (x / chord_length)**4
    )
    
    # Apply a parabolic bending line to the airfoil
    # This creates the "bent" shape
    camber_y = bend_factor * (x / chord_length) * (1 - x / chord_length)
    
    # Combine the thickness and the bending line for upper and lower surfaces
    x_upper = x
    y_upper = y_t + camber_y
    x_lower = x
    y_lower = -y_t + camber_y
    
    # Apply angle of attack rotation
    alpha = math.radians(angle_of_attack)
    cos_alpha = math.cos(alpha)
    sin_alpha = math.sin(alpha)
    
    # Rotate and translate upper surface
    x_upper_rot = x_upper * cos_alpha - y_upper * sin_alpha + x_offset
    y_upper_rot = x_upper * sin_alpha + y_upper * cos_alpha + y_offset
    
    # Rotate and translate lower surface
    x_lower_rot = x_lower * cos_alpha - y_lower * sin_alpha + x_offset
    y_lower_rot = x_lower * sin_alpha + y_lower * cos_alpha + y_offset
    
    return x_upper_rot, y_upper_rot, x_lower_rot, y_lower_rot

def write_airfoils_to_stl(filename, airfoils, span=0.1):
    """
    Generates an STL file with multiple airfoils from a list of parameters.
    """
    with open(filename, 'w') as f:
        f.write("solid three_naca_airfoils\n")
        
        for airfoil_params in airfoils:
            x_upper, y_upper, x_lower, y_lower = generate_bent_airfoil_coordinates(
                n_points=airfoil_params['n_points'],
                chord_length=airfoil_params['chord'],
                x_offset=airfoil_params['x_offset'],
                y_offset=airfoil_params['y_offset'],
                angle_of_attack=airfoil_params['angle'],
                bend_factor=airfoil_params['bend_factor'],
                thickness_front=airfoil_params['thickness_front'],
                thickness_end=airfoil_params['thickness_end']
            )
            
            n = len(x_upper)
            
            # Upper surface triangles
            for i in range(n-1):
                # Triangle 1: (i,0) -> (i+1,0) -> (i,span)
                f.write(f"  facet normal 0.0 0.0 1.0\n")
                f.write(f"    outer loop\n")
                f.write(f"      vertex {x_upper[i]:.6f} {y_upper[i]:.6f} 0.0\n")
                f.write(f"      vertex {x_upper[i+1]:.6f} {y_upper[i+1]:.6f} 0.0\n")
                f.write(f"      vertex {x_upper[i]:.6f} {y_upper[i]:.6f} {span:.6f}\n")
                f.write(f"    endloop\n")
                f.write(f"  endfacet\n")
                
                # Triangle 2: (i+1,0) -> (i+1,span) -> (i,span)
                f.write(f"  facet normal 0.0 0.0 1.0\n")
                f.write(f"    outer loop\n")
                f.write(f"      vertex {x_upper[i+1]:.6f} {y_upper[i+1]:.6f} 0.0\n")
                f.write(f"      vertex {x_upper[i+1]:.6f} {y_upper[i+1]:.6f} {span:.6f}\n")
                f.write(f"      vertex {x_upper[i]:.6f} {y_upper[i]:.6f} {span:.6f}\n")
                f.write(f"    endloop\n")
                f.write(f"  endfacet\n")
            
            # Lower surface triangles
            for i in range(n-1):
                # Triangle 1: (i,0) -> (i,span) -> (i+1,0)
                f.write(f"  facet normal 0.0 0.0 -1.0\n")
                f.write(f"    outer loop\n")
                f.write(f"      vertex {x_lower[i]:.6f} {y_lower[i]:.6f} 0.0\n")
                f.write(f"      vertex {x_lower[i]:.6f} {y_lower[i]:.6f} {span:.6f}\n")
                f.write(f"      vertex {x_lower[i+1]:.6f} {y_lower[i+1]:.6f} 0.0\n")
                f.write(f"    endloop\n")
                f.write(f"  endfacet\n")
                
                # Triangle 2: (i+1,0) -> (i,span) -> (i+1,span)
                f.write(f"  facet normal 0.0 0.0 -1.0\n")
                f.write(f"    outer loop\n")
                f.write(f"      vertex {x_lower[i+1]:.6f} {y_lower[i+1]:.6f} 0.0\n")
                f.write(f"      vertex {x_lower[i]:.6f} {y_lower[i]:.6f} {span:.6f}\n")
                f.write(f"      vertex {x_lower[i+1]:.6f} {y_lower[i+1]:.6f} {span:.6f}\n")
                f.write(f"    endloop\n")
                f.write(f"  endfacet\n")
            
            # End caps (front and back)
            # Front cap (z=0)
            for i in range(n-1):
                if i < n//2: # Upper surface to center
                    f.write(f"  facet normal 0.0 0.0 -1.0\n")
                    f.write(f"    outer loop\n")
                    f.write(f"      vertex {x_upper[i]:.6f} {y_upper[i]:.6f} 0.0\n")
                    f.write(f"      vertex {x_lower[n-1-i]:.6f} {y_lower[n-1-i]:.6f} 0.0\n")
                    f.write(f"      vertex {x_upper[i+1]:.6f} {y_upper[i+1]:.6f} 0.0\n")
                    f.write(f"    endloop\n")
                    f.write(f"  endfacet\n")
            
            # Back cap (z=span)
            for i in range(n-1):
                if i < n//2: # Upper surface to center
                    f.write(f"  facet normal 0.0 0.0 1.0\n")
                    f.write(f"    outer loop\n")
                    f.write(f"      vertex {x_upper[i]:.6f} {y_upper[i]:.6f} {span:.6f}\n")
                    f.write(f"      vertex {x_upper[i+1]:.6f} {y_upper[i+1]:.6f} {span:.6f}\n")
                    f.write(f"      vertex {x_lower[n-1-i]:.6f} {y_lower[n-1-i]:.6f} {span:.6f}\n")
                    f.write(f"    endloop\n")
                    f.write(f"  endfacet\n")
        
        f.write("endsolid three_naca_airfoils\n")

# Parameters for the three airfoils based on the image .
# The 'bend_factor' controls the curvature.
# 'thickness_front' and 'thickness_end' now control the tapering.
airfoils_to_generate = [
    {
        'chord': 2,
        'x_offset': 0.0,
        'y_offset': 0.0,
        'angle': 10.0,
        'n_points': 150,
        'bend_factor': -0.4,
        'thickness_front': 0.25,
        'thickness_end': 0.03
    },
    
    {
        'chord': 0.9,
        'x_offset': 1.75,
        'y_offset': 0.5,
        'angle': 35.0,
        'n_points': 150,
        'bend_factor': -0.3,
        'thickness_front': 0.18,
        'thickness_end': 0.03
    },
    {
        'chord': 0.8,
        'x_offset': 2.35,
        'y_offset': 1.1,
        'angle': 70.0,
        'n_points': 150,
        'bend_factor': -0.2,
        'thickness_front': 0.15,
        'thickness_end': 0.03
    }
]

# Generate the STL file
write_airfoils_to_stl("naca0018.stl", airfoils_to_generate, span=1.0)
print("naca0018.stl file generated successfully with three bent airfoils!")

