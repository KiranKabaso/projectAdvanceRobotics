from sympy import symbols, solve, nsolve
from scipy.spatial.transform import Rotation as Rotation_lib #I want a simple way to get a random rotation matrix

import numpy as np
def find_R_t(points, lines):
  #define k for the line equations c1 + kv1.
  k1, k2, k3 = symbols('k1, k2, k3')
  k_array = [k1, k2, k3]
  # Define R variables
  r11, r12, r13, r21, r22, r23, r31, r32, r33 = symbols('r11 r12 r13 r21 r22 r23 r31 r32 r33')
  r_array = [r11, r12, r13, r21, r22, r23, r31, r32, r33]

  #equations 1 - 9
  #define that Rx + t = c + kV for all points and matching lines for a general 3x3 matrix R.
  equations1_9 = []
  for i in range(3):
    for j in range(i, 3):

      new_equation = ((points[j] - points[i])[0])*r11 + ((points[j] - points[i])[1])*r12 + ((points[j] - points[i])[2])*r13 + (k_array[i]*lines[i][1][0] - k_array[j]*lines[j][1][0]) + (lines[i][0][0] - lines[j][0][0])
      equations1_9.append(new_equation)

      new_equation = ((points[j] - points[i])[0])*r21 + ((points[j] - points[i])[1])*r22 + ((points[j] - points[i])[2])*r23 + (k_array[i]*lines[i][1][1] - k_array[j]*lines[j][1][1]) + (lines[i][0][1] - lines[j][0][1])
      equations1_9.append(new_equation)

      new_equation = ((points[j] - points[i])[0])*r31 + ((points[j] - points[i])[1])*r32 + ((points[j] - points[i])[2])*r33 + (k_array[i]*lines[i][1][2] - k_array[j]*lines[j][1][2]) + (lines[i][0][2] - lines[j][0][2])
      equations1_9.append(new_equation)

  #define that R is a rotation matrix.

  #equations 10-15 -> R orthogonal
  f10 = r11**2 + r12**2 + r13**2 - 1
  f11 = r21**2 + r22**2 + r23**2 - 1
  f12 = r31**2 + r32**2 + r33**2 - 1

  f13 = r11*r21 + r12*r22 + r13*r23
  f14 = r11*r31 + r12*r32 + r13*r33
  f15 = r21*r31 + r22*r32 + r23*r33

  #equation 16 - det(R) = 1
  f16 = (r11*r22*r33 - r11*r23*r32 - r12*r21*r33 + r12*r23*r31 + r13*r21*r32 - r13*r22*r31 - 1)*0

  equations = equations1_9 + [f10, f11, f12, f13, f14, f15, f16]
  variables = k_array + r_array

  #solve
  solution = find_solution(equations, variables)
  if(not solution):
    print("no solution found")
    exit(1)

  #R
  Rrow1 = [solution[3], solution[4], solution[5]]
  Rrow2 = [solution[6], solution[7], solution[8]]
  Rrow3 = [solution[9], solution[10], solution[11]]
  R = [Rrow1, Rrow2, Rrow3]

  #t = c - Rx + kv
  t = lines[0][0] - R @ points[0] + solution[0] * lines[0][1]
  return R, t


def find_solution(equations, variables):
  # solution = solve(equations, variables) #if we want to use solve for a perfect solution, currently failed.
  #using nsolve for a very close solution:
  solution = None
  max_attempts = 1000
  for attempt in range(max_attempts):
      try:
          # Generate a new random initial guess
          initial_guess = np.random.uniform(0, 1, 12)
          # Attempt to solve
          solution = nsolve(equations, variables, initial_guess, tol = 1e-30)
          print(f"Solution found on attempt {attempt + 1}")
          break
      except Exception as e:
          print(f"Attempt {attempt + 1} failed: {e}")
  else:
      print("Failed to find a solution after maximum attempts.")

  return solution


def point_to_line_distance(point, line_point, line_vector):
  diff = point - line_point
  cross_product = np.cross(line_vector, diff)
  distance = np.linalg.norm(cross_product) / np.linalg.norm(line_vector)
  return distance


points = [np.random.rand(3), np.random.rand(3), np.random.rand(3)] #possible solution.

lines_c = [np.random.rand(3), np.random.rand(3), np.random.rand(3)]
lines = [[lines_c[0],points[0] - lines_c[0]], [lines_c[1],points[1] - lines_c[1]], [lines_c[2],points[2] - lines_c[2]]]

rotation_matrix = (Rotation_lib.random().as_matrix())
possible_t = np.random.rand(3)
# possible_t = np.zeros(3)
# rotation_matrix = np.eye(3)
starting_points = [rotation_matrix @ point + possible_t for point in points]

found_R, found_t = find_R_t(starting_points, lines)
if((point_to_line_distance(np.array(found_R@starting_points[0] + found_t, dtype=np.float64), lines[0][0], lines[0][1]) < 1e-15) and (point_to_line_distance(np.array(found_R@starting_points[1] + found_t, dtype=np.float64), lines[1][0], lines[1][1]) < 1e-15) and (point_to_line_distance(np.array(found_R@starting_points[2] + found_t, dtype=np.float64), lines[2][0], lines[2][1])) < 1e-15):
  print('correct')
  print(f"Point Line Distance:  {float(point_to_line_distance(np.array(found_R@starting_points[0] + found_t, dtype=np.float64), lines[0][0], lines[0][1])), float(point_to_line_distance(np.array(found_R@starting_points[1] + found_t, dtype=np.float64), lines[1][0], lines[1][1])), float(point_to_line_distance(np.array(found_R@starting_points[2] + found_t, dtype=np.float64), lines[2][0], lines[2][1]))}")
else: print('incorrect')
# print('original rotation: ',rotation_matrix, 'found rotation: ', found_R)
# print('original translation: ',possible_t, 'found translation: ', found_t)