import numpy as np
import matplotlib.pyplot as plt


# Определение констант
L = 1
alpha_1 = np.pi/6
alpha_3 = 0
c = 3e8
f = 139.466 #139.533
dT_0_squared = 1
T_0_squared = 0.5


# Определение функций
#Лоренцевы контуры
def epsi_11(f):
	# f_j10 = np.array([132.366, 132.433, 134.266, 134.333, 139.466, 139.533])
	f_j10 = np.array([134.266, 134.333, 139.466, 139.533])
	df_j1 = 0.05
	epsi_j10 = np.array([5e-6 / f_j10[i] for i in range(4)])
	return 1 + np.sum((epsi_j10 * (f_j10**2 - f**2) * f_j10 * df_j1) / ((f_j10**2 - f**2)**2 + df_j1**2 * f**2))

def epsi_12(f):
	f_j10 = np.array([132.366, 132.433, 134.266, 134.333, 139.466, 139.533])
	df_j1 = 0.05
	epsi_j10 = np.array([5e-6 / f_j10[i] for i in range(4)])
	return np.sum((epsi_j10 * f_j10**2 * f_j10 * f**2) / ((f_j10**2 - f**2)**2 + df_j1**2 * f**2))

def epsi_31(f):
	f_j30 = np.concatenate((
		132.3 + 0.05 * np.arange(4),
		# 134.2 + 0.05 * np.arange(4),
		# 139.4 + 0.05 * np.arange(4)
	))
	df_j3 = np.array([0.001 * ((i/4) + 3) for i in range(4)])
	epsi_j30 = np.array([0.002 * (5 - (i/4)) * df_j3[i-1] * f_j30[i-1] for i in range(4)])
	asd = 0
	return 1 + np.sum((epsi_j30 * (f_j30**2 - f**2) * f_j30 * df_j3) / ((f_j30**2 - f**2)**2 + df_j3**2 * f**2))

def epsi_32(f):
	f_j30 = np.concatenate((
		132.3 + 0.05 * np.arange(4),
		# 134.2 + 0.05 * np.arange(4),
		# 139.4 + 0.05 * np.arange(4)
	))
	df_j3 = np.array([0.001 * ((i/4) + 3) for i in range(4)])
	epsi_j30 = np.array([0.002 * (5 - (i/4)) * df_j3[i-1] * f_j30[i-1] for i in range(4)])
	return np.sum((epsi_j30 * f_j30**2 * f_j30 * f**2) / ((f_j30**2 - f**2)**2 + df_j3**2 * f**2))

def kz32(f):
	omega = 2 * np.pi * f
	return omega * epsi_32(f) / c * 2 * np.sqrt(epsi_31(f) * np.cos(alpha_1))

def T_sample_squared(f):
	return (T_0_squared + ((np.sqrt(epsi_31(f) * np.cos(alpha_1) / epsi_11(f) * np.cos(alpha_1))) - 1) * dT_0_squared * f) * np.exp(-2 * kz32(f) * L)


# Создание массива значений f и T_sample_squared
f_values = np.arange(139.4, 139.56, 0.05)
T_sample_squared_values = T_sample_squared(f_values)


# Построение графика
plt.plot(f_values, T_sample_squared_values)
plt.xlabel('f')
plt.ylabel('|T_sample|^2')
plt.title('Зависимость |T_sample|^2 от f')
plt.grid(True)
plt.show()