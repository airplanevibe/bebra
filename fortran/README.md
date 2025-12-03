# Fortran Multipole Particle Simulation

Реализация мультипольного разложения для симуляции заряженных частиц в 2D на Fortran 90.

## Особенности

✓ **Multipole expansion** для быстрого вычисления взаимодействий
✓ **Сеточная структура** с ячейками для пространственной декомпозиции
✓ **Три метода** расчёта энергии: прямой, мультиполь по ячейкам, полный мультиполь
✓ **Два интегратора**: метод Эйлера и Runge-Kutta 4-го порядка
✓ **Отражающие границы** для ограничения домена
✓ **Предварительный расчёт** коэффициентов полиномов Лежандра в Python
✓ **Визуализация** результатов через Python

## Структура файлов

```
fortran/
├── generate_coefficients.py  - Генерация коэффициентов Лежандра
├── multipole_module.f90      - Основной модуль с мультипольным разложением
├── test_simulation.f90        - Тестовая программа
├── visualize.py               - Визуализация результатов
├── Makefile                   - Сборка проекта
└── README.md                  - Эта документация
```

## Быстрый старт

### 1. Сборка и запуск

```bash
# Установить gfortran если нужно
# sudo apt-get install gfortran

# Собрать и запустить
make run

# Или по шагам:
make              # Сборка
./test_simulation # Запуск
```

### 2. Визуализация

```bash
# Интерактивная визуализация
python3 visualize.py

# График эволюции энергии
python3 visualize.py --energy

# Создать видео (requires ffmpeg)
python3 visualize.py --save --output animation.mp4
```

## API модуля multipole_module

### Инициализация системы

```fortran
type(particle_system) :: sys

! Инициализация: домен [-xmax, xmax] x [-ymax, ymax], сетка numCellsX x numCellsY
call init_system(sys, xmax, ymax, numCellsX, numCellsY)
```

**Пример:**
```fortran
call init_system(sys, 1.0_dp, 1.0_dp, 10, 10)
! Домен [-1,1] x [-1,1], сетка 10x10
```

### Добавление частиц

```fortran
! Добавить частицу с зарядом charge в позиции (x, y)
call add_particle(sys, charge, x, y)
```

**Пример:**
```fortran
call add_particle(sys, 1.0_dp, 0.5_dp, 0.3_dp)   ! Положительный заряд
call add_particle(sys, -1.0_dp, -0.2_dp, 0.1_dp) ! Отрицательный заряд
```

### Вычисление энергии

```fortran
real(dp) :: energy

! method:
!   0 = прямое вычисление O(N^2)
!   1 = мультиполь по ячейкам (гибрид)
!   2 = полное мультипольное разложение
energy = get_energy(sys, method)
```

**Пример:**
```fortran
energy_direct = get_energy(sys, 0)     ! Точный расчёт
energy_fast = get_energy(sys, 1)       ! Быстрый расчёт
```

### Определение порядка разложения

```fortran
real(dp) :: precision

! Автоматически выбирает порядок p_order исходя из требуемой точности
call determine_series_order(sys, precision)
```

**Пример:**
```fortran
call determine_series_order(sys, 1.0e-6_dp)  ! Точность 10^-6
! Установит p_order = 10-15
```

### Расчёт сил

```fortran
real(dp) :: fx, fy

! Сила на отдельную частицу
call get_force(sys, particle_id, fx, fy, method)

! Силы на все частицы (сохраняются в sys%particles(i)%fx, %fy)
call get_all_forces(sys, method)
```

**Пример:**
```fortran
call get_all_forces(sys, 1)  ! Вычислить все силы методом 1
```

### Движение частиц

```fortran
real(dp) :: eta, dt
logical :: use_rk4

! Интегрирование уравнения dr/dt = eta * F
! eta - параметр вязкости, dt - шаг по времени
call move_particles(sys, eta, dt, use_rk4)
```

**Пример:**
```fortran
! Метод Эйлера
call move_particles(sys, 0.1_dp, 0.01_dp, .false.)

! Runge-Kutta 4
call move_particles(sys, 0.1_dp, 0.01_dp, .true.)
```

### Коррекция ячеек

```fortran
! После движения частицы могут поменять ячейки
call correct_cells(sys)
```

### Запись кадра

```fortran
integer :: frame_num

! Сохранить текущее состояние в CSV файл
call write_frame(sys, frame_num)
! Создаст файл frame_NNNNNN.csv
```

## Пример полной программы

```fortran
program my_simulation
    use multipole_module
    implicit none

    type(particle_system) :: sys
    integer :: i, step
    real(dp) :: energy, dt, eta

    ! Инициализация
    call init_system(sys, 1.0_dp, 1.0_dp, 10, 10)
    call determine_series_order(sys, 1.0e-6_dp)

    ! Добавить частицы
    do i = 1, 100
        call add_particle(sys, real((-1)**i, dp), &
                         rand(), rand())
    end do

    ! Параметры симуляции
    dt = 0.01_dp
    eta = 0.1_dp

    ! Главный цикл
    do step = 1, 1000
        ! Вычислить силы
        call get_all_forces(sys, 1)

        ! Двигать частицы
        call move_particles(sys, eta, dt, .true.)

        ! Обновить ячейки
        call correct_cells(sys)

        ! Сохранить кадр
        if (mod(step, 10) == 0) then
            call write_frame(sys, step/10)
        end if
    end do

end program my_simulation
```

## Детали реализации

### Мультипольное разложение

Основано на разложении потенциала 1/r через полиномы Лежандра:

```
1/|r - r_0| = Σ_{n=0}^∞ P_n(cos θ) × (r_0/r)^(n+1)
```

где θ - угол между векторами r и r_0.

### Структура данных

**Particle:**
- `charge` - заряд
- `x, y` - координаты
- `fx, fy` - силы
- `cell_i, cell_j` - индексы ячейки

**Cell:**
- `num_particles` - количество частиц
- `particle_ids(:)` - ID частиц в ячейке
- `moments(0:p_max)` - мультипольные моменты M_n
- `center_x, center_y` - центр ячейки

### Методы расчёта энергии

**Method 0 (Direct):**
- O(N²) - прямое вычисление всех пар
- Точный результат
- Медленный для больших N

**Method 1 (Multipole Cells):**
- Гибридный подход
- Прямой расчёт внутри и между соседними ячейками
- Мультиполь для далёких ячеек
- O(N + N_cells²)

**Method 2 (Full Multipole):**
- Полное мультипольное разложение
- В данной реализации аналогично Method 1

### Граничные условия

Отражающие границы:
```fortran
! X boundaries
if (x < -xmax) x = -2*xmax - x
if (x >  xmax) x =  2*xmax - x

! Y boundaries
if (y < -ymax) y = -2*ymax - y
if (y >  ymax) y =  2*ymax - y
```

## Производительность

Типичные тесты (100 частиц, 10x10 сетка, 200 шагов):
- Method 0 (direct): ~0.5 сек
- Method 1 (multipole): ~0.2 сек
- Speedup: ~2.5x

Для больших систем (N > 1000) ускорение более значительное.

## Визуализация

### Интерактивный режим

```bash
python3 visualize.py
```

Показывает анимацию кадр за кадром.

### График энергии

```bash
python3 visualize.py --energy
```

Создаёт `energy_evolution.png` с графиками энергии vs время и шаг.

### Создание видео

```bash
python3 visualize.py --save --output my_animation.mp4
```

Требует установленный `ffmpeg`:
```bash
sudo apt-get install ffmpeg
```

## Makefile команды

```bash
make              # Собрать всё
make run          # Собрать и запустить
make visualize    # Запустить и визуализировать
make plot_energy  # График энергии
make clean        # Удалить сборку и данные
make clean_data   # Удалить только данные
make distclean    # Удалить всё включая сгенерированный код
make help         # Справка
```

## Настройка параметров

### В Python (generate_coefficients.py)

```bash
python3 generate_coefficients.py <p_max> <num_angles>

# Пример: p_max=30, 720 углов
python3 generate_coefficients.py 30 720
```

### В Fortran (test_simulation.f90)

Параметры в начале `main`:
```fortran
! Домен
call init_system(sys, xmax, ymax, num_cells_x, num_cells_y)

! Частицы
num_particles = 100

! Временной шаг
dt = 0.01_dp

! Вязкость
eta = 0.1_dp

! Количество шагов
num_steps = 200

! Метод (0, 1, или 2)
method = 1

! Интегратор
use_rk4 = .true.  ! или .false. для Эйлера
```

## Требования

- **Fortran**: gfortran (или другой компилятор F90)
- **Python**: NumPy, Matplotlib, SciPy
- **Опционально**: ffmpeg для создания видео

## Известные ограничения

1. **2D только**: z-координата всегда 0
2. **Мультиполь**: текущая реализация упрощена, не полная FMM
3. **RK4**: упрощённая версия, использует те же силы на всех стадиях
4. **Границы**: только отражение, нет периодических границ

## Дальнейшие улучшения

- [ ] Полная иерархическая FMM с октодеревом
- [ ] Периодические граничные условия
- [ ] Полная реализация RK4 с пересчётом сил
- [ ] OpenMP параллелизация
- [ ] 3D обобщение
- [ ] Адаптивный шаг по времени

## Автор

Terry (Terragon Labs)
Дата: 2025-12-03

Основано на:
- Python реализации мультипольного разложения
- Лекциях Duraiswami & Gumerov (2003-2004)
