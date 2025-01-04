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
  #define S variables
  s11, s12, s13, s21, s22, s23, s31, s32, s33 = symbols('s11 s12 s13 s21 s22 s23 s31 s32 s33')
  s_array = [s11, s12, s13, s21, s22, s23, s31, s32, s33]
  #define P variables
  p1, p2, p3, p4, p5, p6 = symbols('p1 p2 p3 p4 p5 p6')
  p_array = [p1, p2, p3, p4, p5, p6]

  #define r,ij variables
  r2233, r2332, r2133, r2331, r2132, r2231 = symbols('r2233 r2332 r2133 r2331 r2132 r2231')
  rij_array = [r2233, r2332, r2133, r2331, r2132, r2231]
  #define q variables
  q1, q2, q3, q4, q5, q6 = symbols('q1 q2 q3 q4 q5 q6')
  q_array = [q1, q2, q3, q4, q5, q6]
  #define equations
  #1 -43


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

  #equations 10-24 -> R orthogonal
  f10 = r11**2 + r12**2 + r13**2 - 1
  f11 = r21**2 + r22**2 + r23**2 - 1
  f12 = r31**2 + r32**2 + r33**2 - 1

  f13 = s11**2 + s12**2 + s13**2 - 2
  f14 = s21**2 + s22**2 + s23**2 - 2
  f15 = s31**2 + s32**2 + s33**2 - 2

  #S definition equations
  f16 = s11 - r11 - r21
  f17 = s12 - r12 - r22
  f18 = s13 - r13 - r23
  f19 = s21 - r11 - r31
  f20 = s22 - r12 - r32
  f21 = s23 - r13 - r33
  f22 = s31 - r21 - r31
  f23 = s32 - r22 - r32
  f24 = s33 - r23 - r33

  #define det(R) = 1
  #equations 25-43-> det(R) = 1
  f25 = (q1**2 + q2**2 + q3**2 + q4**2 + q5**2 + q6**2) - (p1**2 + p2**2 + p3**2 + p4**2 + p5**2 + p6**2) - 2*(1 + r11**2 + r12**2 + r13**2)

  #p definition equations
  f26 = p1 - 0.5*(r2233**2 -(r22**2 + r33**2))
  f27 = p2 - 0.5*(r2332**2 -(r23**2 + r32**2))
  f28 = p3 - 0.5*(r2133**2 -(r21**2 + r33**2))
  f29 = p4 - 0.5*(r2331**2 -(r23**2 + r31**2))
  f30 = p5 - 0.5*(r2132**2 -(r21**2 + r32**2))
  f31 = p6 - 0.5*(r2231**2 -(r22**2 + r31**2))

  #q definition equations
  f32 = q1 - (r11 + p1)
  f33 = q2 - (r11 - p2)
  f34 = q3 - (r12 - p3)
  f35 = q4 - (r12 + p4)
  f36 = q5 - (r13 + p5)
  f37 = q6 - (r13 - p6)

  #r, ij definition equations
  f38 = r2233 - (r22 + r33)
  f39 = r2332 - (r23 + r32)
  f40 = r2133 - (r21 + r33)
  f41 = r2331 - (r23 + r31)
  f42 = r2132 - (r21 + r32)
  f43 = r2231 - (r22 + r31)

  equations = equations1_9 + [f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22, f23, f24, f25, f26, f27, f28, f29, f30, f31, f32, f33, f34, f35, f36, f37, f38, f39, f40, f41, f42, f43]
  variables = k_array + r_array + s_array + p_array + rij_array + q_array
  #solve the system of equations

  #solve
  solution = find_solution(equations, variables)

  #fill R and find t
  #R
  if(not solution):
    print("no solution found")
    exit(1)
  # if solve
  # Rrow1 = [solution[0][r11], solution[0][r12], solution[0][r13]]
  # Rrow2 = [solution[0][r21], solution[0][r22], solution[0][r23]]
  # Rrow3 = [solution[0][r31], solution[0][r32], solution[0][r33]]
  # R = Matrix([Rrow1, Rrow2, Rrow3])
  # point = Matrix(points[0])
  # #t
  # #t = c - Rx + kv
  # t = lines[0,[0]] - R*point + solution[0][k1] * lines[0,[1]]
  # print(1)
  # print(t)
  # return R, t

  #if nsolve
  Rrow1 = [solution[3], solution[4], solution[5]]
  Rrow2 = [solution[6], solution[7], solution[8]]
  Rrow3 = [solution[9], solution[10], solution[11]]
  R = [Rrow1, Rrow2, Rrow3]
  #t
  #t = c - Rx + kv
  t = lines[0][0] - R @ points[0] + solution[0] * lines[0][1]
  return R, t


def find_solution(equations, variables):
  # solution = solve(equations, variables) #if we want to use solve for a perfect solution, currently failed.

  #use nsolve for approximate solution:
  solution = None
  max_attempts = 100
  for attempt in range(max_attempts):
      try:
          # Generate a new random initial guess
          initial_guess = np.random.uniform(0, 1, 39)
          # Attempt to solve
          solution = nsolve(equations, variables, initial_guess, tol = 1e-11)
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

starting_points = [rotation_matrix @ point + possible_t for point in points]

found_R, found_t = find_R_t(starting_points, lines)
if((point_to_line_distance(np.array(found_R@starting_points[0] + found_t, dtype=np.float64), lines[0][0], lines[0][1]) < 1e-10) and (point_to_line_distance(np.array(found_R@starting_points[1] + found_t, dtype=np.float64), lines[1][0], lines[1][1]) < 1e-10) and (point_to_line_distance(np.array(found_R@starting_points[2] + found_t, dtype=np.float64), lines[2][0], lines[2][1])) < 1e-10):
  print('correct')
  print(f"Point Line Distance:  {float(point_to_line_distance(np.array(found_R@starting_points[0] + found_t, dtype=np.float64), lines[0][0], lines[0][1])), float(point_to_line_distance(np.array(found_R@starting_points[1] + found_t, dtype=np.float64), lines[1][0], lines[1][1])), float(point_to_line_distance(np.array(found_R@starting_points[2] + found_t, dtype=np.float64), lines[2][0], lines[2][1]))}")
else: print('incorrect')
