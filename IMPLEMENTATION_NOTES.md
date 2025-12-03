# Заметки по реализации Grid-Based FMM

## Проблема с текущей реализацией

Текущая "оптимизированная" версия имеет ошибку: она вычисляет far-field вклад в центре ячейки,
но не правильно интерполирует его к частицам.

## Правильный подход для Grid-Based FMM

1. **P2M** (Particle-to-Multipole): Вычислить моменты M_kl для каждой ячейки
2. **M2L** (Multipole-to-Local): ДЛЯ КАЖДОЙ ЧАСТИЦЫ вычислить вклад от далеких ячеек
3. **P2P** (Particle-to-Particle): Прямой расчет для ближних ячеек

Ключевой момент: **НЕ** вычислять локальное разложение в центре ячейки, 
а напрямую вычислять мультипольный вклад для каждой частицы.

## Алгоритм

```
for each particle i:
    phi(i) = 0, fx(i) = 0, fy(i) = 0
    
    # Far field: sum over distant cells
    for each cell j where distance(cell_i, cell_j) > near_range:
        # Evaluate multipole expansion of cell_j at particle i position
        phi(i) += eval_multipole(cell_j.moments, particle_i.position)
        f(i) += grad_multipole(cell_j.moments, particle_i.position)
    
    # Near field: direct sum over nearby particles
    for each cell j where distance(cell_i, cell_j) <= near_range:
        for each particle k in cell_j:
            if i != k:
                phi(i) += 1/|r_i - r_k|
                f(i) += (r_i - r_k)/|r_i - r_k|^3
```

Это Grid-Based FMM без M2L в классическом смысле - мы не делаем трансляцию мультиполей,
а напрямую вычисляем их в точках частиц.

## Преимущества

- Простота реализации
- Точность (никаких аппроксимаций при трансляции)
- Все еще O(N) для N >> ngrid^2

## Сложность

- P2M: O(N * p^2)
- Far-field: O(N * ngrid^2 * p^2)
- Near-field: O(N * particles_per_cell * near_range^2)

Для ngrid = O(N^{1/2}) и near_range = O(1):
- Total: O(N * N * p^2) = O(N^{3/2} * p^2)

Лучше чем O(N^2) но не идеальный O(N).

Для истинного O(N) нужна полная M2L трансляция с Taylor expansion.

