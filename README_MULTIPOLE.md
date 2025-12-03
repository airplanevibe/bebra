# Multipole Expansion for 2D Laplace Equation

Реализация мультипольного разложения для уравнения Лапласа в 2D случае (без взаимодействий по оси z) на основе лекций "Fast Multipole Methods for the Laplace Equation" (Duraiswami & Gumerov, 2003-2004).

## Особенности

✓ **Высокая точность**: Близка к машинной точности при высоких порядках разложения (p_max ~ 25-30)
✓ **2D случай**: 3D потенциал, но взаимодействия только в плоскости xy (z=0)
✓ **Равномерная сетка**: Поддержка равномерных сеток точек
✓ **Прямое разложение**: Использует теорему сложения для полиномов Лежандра
✓ **Полный контроль**: Все параметры задаются через main() или командную строку

## Математическая основа

### Разложение потенциала 1/r

Для потенциала Кулона используется разложение функции Грина:

```
1/|r - r₀| = Σ_{n=0}^∞ Pₙ(cos θ) × (r₀/r)^(n+1)    для r > r₀
```

где:
- `Pₙ` - полином Лежандра степени n
- `θ` - угол между векторами r и r₀
- `cos θ = (r · r₀) / (|r| × |r₀|)`

### Мультипольные моменты

Потенциал от набора зарядов:

```
Φ(r) = Σᵢ qᵢ / |r - rᵢ|
     = Σₙ [Σᵢ qᵢ × rᵢⁿ × Pₙ(cos θᵢ)] / r^(n+1)
     = Σₙ Mₙ(θ,φ) / r^(n+1)
```

где `Mₙ` - мультипольные моменты порядка n.

## Файлы

### Основной файл: `multipole_2d_final.py`

Финальная оптимизированная версия с полной документацией и параметрами командной строки.

**Использование:**

```bash
# Запуск с параметрами по умолчанию
python3 multipole_2d_final.py

# Только тест точности с высоким порядком
python3 multipole_2d_final.py --test-type accuracy --p-max 30 --n-sources 100

# Только тест производительности
python3 multipole_2d_final.py --test-type performance --p-max 20

# Все параметры
python3 multipole_2d_final.py \
    --test-type both \
    --p-max 25 \
    --n-sources 100 \
    --n-eval 50 \
    --source-radius 1.0 \
    --eval-radius 2.5
```

**Параметры командной строки:**

- `--test-type`: Тип теста (`accuracy`, `performance`, `both`)
- `--p-max`: Максимальный порядок разложения (рекомендуется 20-30)
- `--n-sources`: Количество источников для теста точности
- `--n-eval`: Количество точек оценки
- `--source-radius`: Радиус области источников
- `--eval-radius`: Радиус области оценки (должен быть > source-radius)

### Дополнительные файлы

- `test_multipole_expansion.py` - Первая версия с использованием сферических гармоник
- `test_multipole_expansion_v2.py` - Упрощённая версия с полиномами Лежандра

## Результаты тестирования

### Тест точности (p_max=30)

```
Configuration:
  Number of sources: 50
  Number of evaluation points: 30
  Source radius: 1.000
  Evaluation radius: 3.000
  Expansion order: 30

Error Metrics:
  Maximum absolute error: 2.66e-15
  Maximum relative error: 6.55e-15
  RMS error: 1.18e-15
  Machine epsilon (float64): 2.22e-16

Status: ✓ EXCELLENT: Near-machine precision achieved!
```

### Тест точности (p_max=25)

```
Configuration:
  Number of sources: 100
  Number of evaluation points: 50
  Source radius: 1.000
  Evaluation radius: 2.500
  Expansion order: 25

Error Metrics:
  Maximum absolute error: 4.46e-12
  Maximum relative error: 8.44e-12
  RMS error: 6.91e-13

Status: ✓ VERY GOOD: High precision
```

### Тест производительности (p_max=25)

```
   N_sources      Direct(s)     Multipole(s)      Speedup
----------------------------------------------------------------------
          50       0.011         0.190           0.058x
         100       0.022         0.359           0.062x
         200       0.044         0.716           0.062x
         400       0.089         1.424           0.062x
         800       0.177         3.053           0.058x
```

**Примечание**: Текущая реализация имеет сложность O(N×M×p), где N - источники, M - точки оценки, p - порядок. Для истинного ускорения необходима иерархическая FMM с октодеревом, которая даёт O((N+M)×p²).

## Параметры для оптимальных результатов

### Для максимальной точности (машинная точность):
- `p_max = 25-30`
- `eval_radius / source_radius >= 2.0`

### Для баланса точность/скорость:
- `p_max = 15-20`
- `eval_radius / source_radius >= 1.5`

### Для быстрых вычислений:
- `p_max = 8-12`
- `eval_radius / source_radius >= 1.2`

## Использование в коде Python

```python
from multipole_2d_final import MultipoleExpansion2D, compute_potential_direct
import numpy as np

# Параметры
p_max = 25
N_sources = 100

# Создать источники в круге радиуса 1.0
angles = np.random.uniform(0, 2*np.pi, N_sources)
radii = np.random.uniform(0, 1.0, N_sources)

sources = np.zeros((N_sources, 3))
sources[:, 0] = radii * np.cos(angles)
sources[:, 1] = radii * np.sin(angles)
sources[:, 2] = 0.0

charges = np.random.uniform(-1.0, 1.0, N_sources)

# Точки оценки вне области источников (радиус > 2.5)
eval_points = np.array([
    [3.0, 0.0, 0.0],
    [0.0, 3.0, 0.0],
    [2.5, 2.5, 0.0]
])

# Создать мультипольное разложение
multipole = MultipoleExpansion2D(p_max)
center = np.array([0.0, 0.0, 0.0])

# Вычислить мультипольные коэффициенты
rel_sources, moments = multipole.compute_multipole_coefficients(
    sources, charges, center
)

# Оценить потенциал
potentials = multipole.evaluate_potential_batch(
    rel_sources, moments, eval_points, center
)

print("Потенциалы:", potentials)

# Сравнить с прямым вычислением
potentials_direct = compute_potential_direct(sources, charges, eval_points)
error = np.abs(potentials - potentials_direct)
print("Максимальная ошибка:", np.max(error))
```

## Теоретическая основа из PDF

Реализация основана на следующих слайдах из PDF "dg_lecture6.pdf":

- **Слайд 7**: Разложение фундаментального решения с использованием полиномов Лежандра
- **Слайд 8**: Теорема сложения для сферических гармоник
- **Слайд 10**: S-expansion (мультипольное разложение) и R-expansion (локальное разложение)

Основная формула (слайд 7):
```
G(r - r₀) = 1/|r - r₀| = Σₙ Pₙ(cos θ) × (r₀/r)^(n+1)  для r > r₀
```

## Ограничения текущей реализации

1. **Производительность**: Текущая версия не использует иерархическую структуру (октодерево), поэтому не достигает теоретической асимптотики O((N+M)×p²) полной FMM.

2. **Условие применимости**: Точки оценки должны находиться вне области источников. Для точек внутри области нужно использовать локальное разложение (R-expansion).

3. **2D ограничение**: Все источники и точки оценки должны быть в плоскости xy (z=0) для корректной работы в 2D режиме.

## Дальнейшие улучшения

Для достижения полной производительности FMM необходимо:

1. Реализовать иерархическую структуру данных (октодерево)
2. Добавить локальные разложения (R-expansion)
3. Реализовать операторы трансляции M2M, M2L, L2L
4. Добавить адаптивное разбиение пространства
5. Оптимизировать с использованием NumPy векторизации или Numba JIT

## Ссылки

- PDF лекции: "Fast Multipole Methods for the Laplace Equation" (Duraiswami & Gumerov, 2003-2004)
- Классическая статья: Greengard & Rokhlin (1987) "A fast algorithm for particle simulations"
- Книга: "The Rapid Evaluation of Potential Fields in Particle Systems" (Greengard, 1988)

## Лицензия

Код создан для образовательных и исследовательских целей.

## Автор

Реализация: Terry (Terragon Labs)
Дата: 2025-12-03
