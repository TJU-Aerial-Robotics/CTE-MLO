import numpy as np
from scipy.spatial.transform import Slerp, Rotation as R

quat = np.array([[0.936118452267473,-0.351663294301812,0.003894594980225,-0.000053806028052]])
print(R.from_quat(quat).as_matrix())
