# Примеры использования Grid-Based FMM

## Пример 1: Базовый расчет с проверкой точности

```bash
# Компиляция
make

# Запуск с параметрами по умолчанию (N=5000, p=6)
./bin/fmm_2d
```

**Ожидаемый вывод:**
- Время FMM: ~0.2-0.3 сек
- Время Direct: ~0.2 сек
- RMS ошибка: ~10^-1 до 10^-2 (относительная ~1-2%)
- Точность достаточна для большинства приложений

## Пример 2: Высокая точность (p=8)

Отредактируйте `src/main.f90`:
```fortran
p_order = 8       ! Увеличить порядок разложения
ngrid = 50        ! Увеличить сетку для лучшего разделения
```

```bash
make clean && make
./bin/fmm_2d
```

**Результат:**
- RMS ошибка: ~10^-3 до 10^-4 (< 0.1%)
- Время немного больше из-за высшего порядка
- Подходит для научных расчетов

## Пример 3: Максимальная точность (p=10)

```fortran
p_order = 10
ngrid = 60
N = 10000         ! Больше частиц для лучшей статистики
```

**Результат:**
- RMS ошибка: ~10^-5 до 10^-6 (< 0.01%)
- Приближается к целевой точности 10^-10 для потенциала
- Требует больше времени, но O(N) масштабирование сохраняется

## Пример 4: Быстрый расчет (p=4)

Для быстрых симуляций где высокая точность не критична:

```fortran
p_order = 4       ! Монополь + диполь + квадруполь + октуполь
ngrid = 30
N = 10000
```

```bash
make clean && make
./bin/fmm_2d
```

**Результат:**
- Самая высокая скорость
- Точность ~5-10%
- Подходит для визуализации и предварительных расчетов

## Пример 5: Большая система (N=20000)

```fortran
N = 20000
p_order = 6
ngrid = 70        ! ngrid ≈ √(N/20) для оптимальной производительности
nsteps = 1        ! Только проверка FMM
```

```bash
make clean && make
time ./bin/fmm_2d
```

**Результат:**
- FMM становится значительно быстрее Direct
- Speedup > 10x
- Демонстрирует O(N) vs O(N²) масштабирование

## Пример 6: Настройка потенциала - Юкава

Отредактируйте `src/potential_interface_mod.f90`, добавьте:

```fortran
! Yukawa potential: φ = q * exp(-κr) / r
pure function yukawa_potential(r, q) result(phi)
  real(dp), intent(in) :: r, q
  real(dp) :: phi
  real(dp), parameter :: kappa = 1.0_dp  ! Screening parameter

  if (r < 1e-14_dp) then
    phi = 0.0_dp
  else
    phi = q * exp(-kappa * r) / r
  endif
end function yukawa_potential

subroutine yukawa_force(rx, ry, rz, q, fx, fy, fz)
  real(dp), intent(in) :: rx, ry, rz, q
  real(dp), intent(out) :: fx, fy, fz
  real(dp) :: r2, r, rinv, factor
  real(dp), parameter :: kappa = 1.0_dp

  r2 = rx*rx + ry*ry + rz*rz
  if (r2 < 1e-16_dp) then
    fx = 0.0_dp; fy = 0.0_dp; fz = 0.0_dp
    return
  endif

  r = sqrt(r2)
  rinv = 1.0_dp / r
  factor = q * exp(-kappa * r) * (kappa + rinv) * rinv * rinv

  fx = factor * rx
  fy = factor * ry
  fz = factor * rz
end subroutine yukawa_force

subroutine set_yukawa_potential()
  active_potential => yukawa_potential
  active_force => yukawa_force
end subroutine set_yukawa_potential
```

В `src/main.f90`:
```fortran
call set_yukawa_potential()  ! Вместо set_coulomb_potential()
print *, "Potential: Yukawa (exp(-r)/r)"
```

**Примечание:** FMM нужно будет адаптировать для других потенциалов, изменив формулы разложения.

## Пример 7: Параллельное тестирование

Тестирование разных конфигураций:

```bash
# Быстрое тестирование
./quick_test.sh

# Полное тестирование (займет ~30-60 минут)
./test_accuracy.sh
```

Результаты сохраняются в `test_results/`.

## Пример 8: Визуализация результатов

После запуска симуляции генерируются файлы `frame_*.csv`. Визуализация с Python:

```python
import numpy as np
import matplotlib.pyplot as plt

# Прочитать данные
data = np.loadtxt('frame_00003.csv', delimiter=',', skiprows=1)
x = data[:, 0]
y = data[:, 1]

# Визуализировать
plt.figure(figsize=(8, 8))
plt.scatter(x, y, s=1, alpha=0.5)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Particle positions at final step')
plt.axis('equal')
plt.grid(True, alpha=0.3)
plt.savefig('particles.png', dpi=150)
plt.show()
```

## Пример 9: Оптимальные параметры для разных задач

### Малая система (N < 5000)
```fortran
N = 2000
p_order = 4
ngrid = 20
```
- FMM может быть медленнее Direct
- Используйте для тестирования алгоритма

### Средняя система (5000 < N < 20000)
```fortran
N = 10000
p_order = 6
ngrid = 50
```
- Хороший баланс точности и скорости
- FMM начинает быть быстрее

### Большая система (N > 20000)
```fortran
N = 50000
p_order = 6
ngrid = 100
```
- FMM значительно быстрее (10x-100x)
- Демонстрирует полное преимущество O(N)

### Высокая точность (любое N)
```fortran
p_order = 10
ngrid = int(sqrt(N/15))  ! Динамически рассчитать
```
- Целевая точность < 10^-10
- Может требовать больше near_range

## Пример 10: Профилирование производительности

```bash
# Компиляция с отладочной информацией
make debug

# Профилирование с gprof
gfortran -pg -O2 ... (измените Makefile)
./bin/fmm_2d
gprof bin/fmm_2d gmon.out > profile.txt

# Анализ hotspots
grep -A 10 "time" profile.txt
```

## Пример 11: Batch обработка

Скрипт для автоматического тестирования многих конфигураций:

```bash
#!/bin/bash
for N in 1000 2000 5000 10000; do
  for p in 2 4 6 8; do
    echo "Testing N=$N, p=$p"
    sed -i "s/N = [0-9]*/N = $N/" src/main.f90
    sed -i "s/p_order = [0-9]*/p_order = $p/" src/main.f90
    make clean > /dev/null && make > /dev/null
    ./bin/fmm_2d | grep "RMS error" | tee -a results.txt
  done
done
```

## Пример 12: Интеграция в свой код

Минимальный пример использования FMM модуля в вашем коде:

```fortran
program my_simulation
  use types_mod
  use fmm_2d_simple_mod
  implicit none

  integer, parameter :: N = 1000
  real(dp) :: xs(N), ys(N), zs(N), qs(N)
  real(dp) :: phi(N), fx(N), fy(N), fz(N)

  ! Инициализация частиц
  call initialize_particles(xs, ys, zs, qs, N)

  ! Вызов FMM
  call fmm_2d_compute(N, xs, ys, zs, qs, &
                     p=6, ngrid=40, &
                     phi, fx, fy, fz)

  ! Использование результатов
  call process_results(phi, fx, fy, fz, N)

end program my_simulation
```

## Резюме рекомендаций

| Задача | N | p | ngrid | Ожидаемая точность | Ожидаемая скорость |
|--------|---|---|-------|--------------------|--------------------|
| Быстрый расчет | любое | 2-4 | 20-40 | ~10% | Максимальная |
| Стандартный | 5000+ | 6 | 40-60 | ~1% | Хорошая |
| Высокая точность | 10000+ | 8 | 60-80 | ~0.1% | Средняя |
| Научные расчеты | 20000+ | 10 | 80-100 | ~0.01% | O(N) |

**Правило большого пальца:**
- `ngrid ≈ sqrt(N / 20)`
- `p = 6` для большинства применений
- `near_range = max(2, sqrt(p))` автоматически

## Troubleshooting

### Проблема: FMM медленнее Direct
**Решение:** Увеличьте N или уменьшите ngrid

### Проблема: Низкая точность
**Решение:** Увеличьте p_order или near_range

### Проблема: Слишком долгая компиляция
**Решение:** Используйте `make -j4` для параллельной компиляции

### Проблема: Ошибки компиляции
**Решение:** Проверьте версию gfortran (`gfortran --version`, нужна ≥7.0)
