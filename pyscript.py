import os
import imageio.v2 as imageio
import numpy as np
import matplotlib.pyplot as plt
import glob

# Сортируем файлы по порядку
filenames = sorted(glob.glob("bin/animation/velocity_out_*.bin"))

# Создаём пустой список для изображений
images = []

# Минимальное и максимальное значения p для фиксирования оси
p_min = float('inf')
p_max = float('-inf')

# Перебираем файлы дважды - первый раз для нахождения
# максимальных и минимальных значений p
for filename in filenames:
    # Читаем бинарный файл с использованием numpy
    data = np.fromfile(filename, dtype=np.float32)
    # Переформатируем данные в двухмерный массив (т.к. у вас две колонки)
    data = data.reshape(-1, 2)
    p_min = min(p_min, np.min(data[:, 1]))
    p_max = max(p_max, np.max(data[:, 1]))

for filename in filenames:
    # Читаем бинарный файл с использованием numpy
    data = np.fromfile(filename, dtype=np.float32)

    # Переформатируем данные в двухмерный массив (т.к. у вас две колонки)
    data = data.reshape(-1, 2)

    # Строим график
    plt.figure(figsize=(10, 5))
    plt.plot(data[:, 0], data[:, 1])
    plt.title(filename)
    plt.xlabel("x")
    plt.ylabel("p")
    plt.ylim([p_min, p_max])
    plt.grid(True)

    # Сохраняем график в файл
    image_filename = f"{filename}.png"
    plt.savefig(image_filename)

    # Закрываем график
    plt.close()

    # Добавляем изображение в список
    images.append(imageio.imread(image_filename))

    # Удаляем временный файл изображения
    os.remove(image_filename)

# Создаём gif из изображений
# Включаем зацикливание с помощью loop=0
imageio.mimsave('output.gif', images, duration=0.5, loop=0)
