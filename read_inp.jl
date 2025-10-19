### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 9bc01b20-8bc2-11f0-32ac-a587080cee1e
begin

	####################### ЧИТКА ИЗ ФАЙЛА ##############################
	using DelimitedFiles
	using LinearAlgebra
	using Statistics 
	
	function parse_abaqus_file(filename)
	    nodes_list = []  # Список для хранения узлов
	    elements_list = []  # Список для хранения элементов
	    reading_nodes = false
	    reading_elements = false
	
	    open(filename, "r") do file
	        for line in eachline(file)
	            line = strip(line)
	            if isempty(line)
	                continue
	            end
	
	            if startswith(line, "*NODE")
	                reading_nodes = true
	                reading_elements = false
	                continue
	            elseif startswith(line, "*ELEMENT")
	                reading_elements = true
	                reading_nodes = false
	                continue
	            end
	
	            if reading_nodes
	                parts = split(line, ',')
	                if length(parts) < 3
	                    continue
	                end

	                x = parse(Float64, parts[2])
	                y = parse(Float64, parts[3])
	                push!(nodes_list, (x, y))
	            elseif reading_elements
	                parts = split(line, ',')
	                if length(parts) < 5
	                    continue
	                end
	                n1 = parse(Int, parts[2])
	                n2 = parse(Int, parts[3])
	                n3 = parse(Int, parts[4])
	                n4 = parse(Int, parts[5])
	                push!(elements_list, [n1, n2, n3, n4])
	            end
	        end
	    end

	    nodes_matrix = zeros(length(nodes_list), 2)
	    for (i, (x, y)) in enumerate(nodes_list)
	        nodes_matrix[i, :] = [x, y]
	    end

	    connectivity_matrix = zeros(Int, length(elements_list), 4)
	    for (i, conn) in enumerate(elements_list)
	        connectivity_matrix[i, :] = conn
	    end
	
	    xs = nodes_matrix[:, 1]
	    ys = nodes_matrix[:, 2]
	    unique_x = sort(unique(round.(xs, digits=10)))
	    unique_y = sort(unique(round.(ys, digits=10)))
	    Lx = maximum(unique_x)
	    Ly = maximum(unique_y)
	    nx = length(unique_x) - 1
	    ny = length(unique_y) - 1
	
	    return Lx, Ly, nx, ny, nodes_matrix, connectivity_matrix
	end
	

	Lx, Ly, nx, ny, nodes, connectivity = parse_abaqus_file("C:\\Users\\Admin\\res2.inp")
	println("Lx = ", Lx)
	println("Ly = ", Ly)
	println("nx = ", nx)
	println("ny = ", ny)

	println("\nNodes (координаты x, y):")
	for (i, row) in enumerate(eachrow(nodes))
    	println("Узел $i: ($(row[1]), $(row[2]))")
	end

	println("\nConnectivity (элементы и их узлы):")
	println("Элемент\tУзел1\tУзел2\tУзел3\tУзел4")
	for i in 1:size(connectivity, 1)
    	println("$i\t$(connectivity[i, 1])\t$(connectivity[i, 2])\t$(connectivity[i, 3])\t$(connectivity[i, 4])")
	end

	##################################################################################################################################################################
	# ********************** ПАРАМЕТРЫ ЗАДАЧИ **********************
E = 70000 # Модуль Юнга (МПа) - характеристика жесткости материала
nu = 0.3 # Коэффициент Пуассона - характеризует поперечную деформацию
thickness = 0.1 # Толщина пластины (м) 
pressure = 100 # Давление (МПа) - нагрузка на правую грань пластины

#Lx = 100 # Длина пластины по оси x (мм)
#Ly = 20 # Длина пластины по оси y (мм)

#nx = 100 #8   #20   #50     # Число элементов по длине (оси x)
#ny = 20 #1   #4    #10     # Число элементов по ширине (оси y)

# Расчет общего количества узлов в сетке
num_nodes = (nx + 1) * (ny + 1) # Формула: (элементы+1) по x * (элементы+1) по y

# ********************** МАТЕРИАЛЬНЫЕ ХАРАКТЕРИСТИКИ **********************

# Функция возвращает матрицу упругости D для плоского напряженного состояния
function get_D_matrix()
    # матрица упругости для изотропного материала: ПНС
    D1 = E / (1 - nu * nu) * [
        1     nu     0;
        nu     1     0;
        0     0   (1-nu)/2
    ]
    factor = E / ((1 + nu) * (1 - 2 * nu))
    D2 = factor * [
        1 - nu     nu         0;
        nu       1 - nu       0;
        0          0     (1 - 2 * nu)/2
    ]
    return D2
end

# ********************** ФУНКЦИИ ФОРМЫ И ИХ ПРОИЗВОДНЫЕ **********************

# Функции формы для 4-узлового изопараметрического элемента
# xi, eta - естественные координаты в диапазоне [-1, 1]
function shape_functions(xi, eta) 
    # Вектор из 4 функций формы
	N = 0.25 * [
        (1 - xi)*(1 - eta), # Для первого узла
        (1 + xi)*(1 - eta), # Для второго узла
        (1 + xi)*(1 + eta), # Для третьего узла
        (1 - xi)*(1 + eta)  # Для четвертого узла
    ]
    return N
end

# Производные функций формы по естественным координатам xi и eta
function shape_function_derivatives(xi, eta) 
    # Матрица 4x2 (4 узла, 2 производные для каждого)
    dNdxi = zeros(4, 2)

    # Производные для первого узла
    dNdxi[1, 1] = -0.25 * (1 - eta) # ∂N₁/∂ξ
    dNdxi[1, 2] = -0.25 * (1 - xi)  # ∂N₁/∂η

    # Производные для второго узла
    dNdxi[2, 1] = 0.25 * (1 - eta)  # ∂N₂/∂ξ
    dNdxi[2, 2] = -0.25 * (1 + xi) # ∂N₂/∂η

    # Производные для третьего узла
    dNdxi[3, 1] = 0.25 * (1 + eta)  # ∂N₃/∂ξ
    dNdxi[3, 2] = 0.25 * (1 + xi)  # ∂N₃/∂η

    # Производные для четвертого узла
    dNdxi[4, 1] = -0.25 * (1 + eta) # ∂N₄/∂ξ
    dNdxi[4, 2] = 0.25 * (1 - xi)   # ∂N₄/∂η

    return dNdxi
end

# ********************** ГЕОМЕТРИЧЕСКИЕ ПРЕОБРАЗОВАНИЯ **********************/

# Вычисление якобиана преобразования между естественными и глобальными координатами
function jacobian(dNdxi, node_coords) 
    J = zeros(2, 2) # Инициализация нулевой матрицы 2x2

    # Суммирование вкладов от всех узлов элемента:
	for i in 1:4
        J[1,1] += dNdxi[i, 1] * node_coords[i, 1] # ∂x/∂ξ
        J[1,2] += dNdxi[i, 1] * node_coords[i, 2] # ∂y/∂ξ
        J[2,1] += dNdxi[i, 2] * node_coords[i, 1] # ∂x/∂η
        J[2,2] += dNdxi[i, 2] * node_coords[i, 2] # ∂y/∂η
    end
    return J
end

# ********************** МАТРИЦА ДЕФОРМАЦИЙ **********************/

# Формирование матрицы B (связь деформаций с перемещениями) для одного узла
function B_matrix(dNdx, node)
	# Матрица 3x2 для плоской задачи
	B = [
        dNdx[node, 1] 0;
        0 dNdx[node, 2];
        dNdx[node, 2] dNdx[node, 1]
    ]
    return B
end

# ********************** СБОРКА ГЛОБАЛЬНОЙ СИСТЕМЫ **********************/

# Сборка глобальной матрицы жесткости из матриц элементов
function assemble_global_stiffness(element_stiffness, connectivity)
    # Инициализация глобальной матрицы жесткости
    # Размер: 2*число_узлов x 2*число_узлов (по 2 степени свободы на узел)
    K = zeros(2 * num_nodes, 2 * num_nodes)

    # Цикл по всем элементам
	for e in 1:size(element_stiffness, 1)
    # Цикл по узлам элемента (i - локальный номер узла в элементе)
        for i in 1:4
            ni = connectivity[e, i] # Глобальный номер i-го узла

            # Цикл по узлам элемента (j - локальный номер узла в элементе)
            for j in 1:4
                nj = connectivity[e, j] # Глобальный номер j-го узла

                # Добавление вклада от элемента в глобальную матрицу
                # выделяет блок 2x2 (по 2 степени свободы на узел)
				K[2*(ni-1)+1 : 2*(ni-1)+2, 2*(nj-1)+1 : 2*(nj-1)+2] .+= element_stiffness[e][2*i-1 : 2*i, 2*j-1 : 2*j]
			end
		end
	end
    return K
end

# Функция для вычисления напряжений в элементах
function calculate_stresses(element_stiffness, connectivity, U, nodes) 
    # Вектор для хранения напряжений в каждом элементе (в виде векторов напряжений)
    element_stresses = [zeros(3) for _ in 1:size(connectivity, 1)]
    
    # Матрица упругости (уже вычислена ранее)
    D = get_D_matrix()
    
    # Цикл по всем элементам для вычисления напряжений
    for e in 1:size(connectivity, 1)
        # 1. Получаем координаты узлов текущего элемента
        node_coords = zeros(4, 2)
        for i in 1:4
			node_coords[i, 1] = nodes[connectivity[e, i], 1] # x-координата
			node_coords[i, 2] = nodes[connectivity[e, i], 2]  # y-координата
		end
        
        # 2. Получаем перемещения узлов текущего элемента
        element_displacements = zeros(8) # 4 узла * 2 степени свободы
        for i in 1:4
            node = connectivity[e, i]
            element_displacements[2*i-1] = U[2*node-1] # x-перемещение
            element_displacements[2*i] = U[2*node] # y-перемещение
        end
        
        # 3. Вычисляем напряжения в центре элемента (xi=0, eta=0)
        xi = 0.0
		eta = 0.0
        
        # 3.1. Вычисляем производные функций формы
        dNdxi = shape_function_derivatives(xi, eta)
        
        # 3.2. Вычисляем якобиан
        J = jacobian(dNdxi, node_coords)
        detJ = det(J)
        
        # 3.3. Вычисляем производные по глобальным координатам
        dNdx=zeros(4, 2)
        invJ = inv(J)
        
        for i in 1:4
            dNdx[i, 1] = invJ[1, 1]*dNdxi[i, 1] + invJ[1, 2]*dNdxi[i, 2] # dN/dx
            dNdx[i, 2] = invJ[2, 1]*dNdxi[i, 1] + invJ[2, 2]*dNdxi[i, 2] # dN/dy
		end
        
        # 3.4. Формируем матрицу B для элемента
        B = zeros(3, 8);
        for i in 1:4
            Bi = B_matrix(dNdx, i)
			
			B[:, 2*i - 1 : 2*i] .= Bi
        end
        
        # 3.5. Вычисляем деформации: ε = B * u
        strains = B * element_displacements
        
        # 3.6. Вычисляем напряжения: σ = D * ε
        stresses = D * strains
        
        # 4. Сохраняем напряжения для текущего элемента
        element_stresses[e] = stresses
	end
    
    return element_stresses
end

	
# РАСКОММЕНТИРОВАТЬ, ЕСЛИ ХОТИМ ИСПОЛЬЗОВАТЬ СВОЮ СЕТКУ, А НЕ ЧИТАТЬ ЕЕ ИЗ ФАЙЛА
	
# ***** 1. ГЕНЕРАЦИЯ КОНЕЧНОЭЛЕМЕНТНОЙ СЕТКИ *****/
# Координаты всех узлов сетки (num_nodes x 2)
 num_nodes = (nx + 1) * (ny + 1)
#nodes = zeros(num_nodes, 2)  # Массив размером num_nodes × 2 для хранения координат узлов
# Равномерная сетка — узлы располагаются на равных расстояниях
#for j in 0:ny
#    for i in 0:nx
#        node_id = j * (nx + 1) + i + 1  # индексация с 1
#        nodes[node_id, 1] = Lx * i / nx  # x-координата
#        nodes[node_id, 2] = Ly * j / ny  # y-координата
#    end
#end

#***** 2. СОЗДАНИЕ ТАБЛИЦЫ СВЯЗНОСТИ *****/

# Каждый элемент связывает 4 узла (квадратные элементы)
#	function generate_connectivity(nx, ny)
#    conn = zeros(Int, nx * ny, 4)
            # Нумерация узлов элемента против часовой стрелки:
            # 3 ---- 2
            # |      |
            # |      |
            # 0 ---- 1
#    for j in 0:(ny - 1)
#        for i in 0:(nx - 1)
#            elem_id = j * nx + i + 1
#            conn[elem_id, :] = [
#                j * (nx + 1) + i + 1,
#                j * (nx + 1) + i + 2,
#                (j + 1) * (nx + 1) + i + 2,
#                (j + 1) * (nx + 1) + i + 1
#            ]
#        end
#    end

#    return conn
#end
#	connectivity = generate_connectivity(nx, ny)


	# ***** 3. МАТРИЦА УПРУГОСТИ D *****
    D = get_D_matrix()

    # ***** 4. ПАРАМЕТРЫ ИНТЕГРИРОВАНИЯ ПО ГАУССУ *****

    # Точки интегрирования (2 точки Гаусса в каждом направлении)
     gauss_points = [ -1 / sqrt(3) 1 / sqrt(3) ]
	# Весовые коэффициенты (для 2 точек равны 1)
     gauss_weights = [ 1.0 1.0 ]

    #gauss_points = [ -sqrt(3.0 / 5.0) 0 sqrt(3.0 / 5.0) ]
    #gauss_weights = [ 5.0 / 9.0 8.0 / 9.0 5.0 / 9.0 ]

	# ***** 5. ВЫЧИСЛЕНИЕ МАТРИЦ ЖЕСТКОСТИ ЭЛЕМЕНТОВ *****
# Инициализируем массив матриц жесткости для всех элементов (по 8 степеней свободы)
element_stiffness = [zeros(8, 8) for _ in 1:size(connectivity, 1)]

# Цикл по всем элементам
for e in 1:size(connectivity, 1)
    # Получаем координаты узлов текущего элемента
    node_coords = zeros(4, 2)
    for i in 1:4
        node_id = connectivity[e, i]
        node_coords[i, 1] = nodes[node_id, 1] # x-координата
        node_coords[i, 2] = nodes[node_id, 2] # y-координата
    end

    # Интегрирование по площади элемента (Гаусс 3 точки)
    for gp_xi in eachindex(gauss_points), gp_eta in eachindex(gauss_weights)
        xi = gauss_points[gp_xi]
        eta = gauss_points[gp_eta]
        weight = gauss_weights[gp_xi] * gauss_weights[gp_eta]

        # 5.1 Вычисление производных функций формы
        dNdxi = shape_function_derivatives(xi, eta)

        # 5.2 Вычисление якобиана
        J = jacobian(dNdxi, node_coords)
        detJ = det(J)

        # 5.3 Преобразование производных в глобальные координаты
        invJ = inv(J)
        dNdx = zeros(4, 2)
        for i in 1:4
            dNdx[i, 1] = invJ[1, 1]*dNdxi[i, 1] + invJ[1, 2]*dNdxi[i, 2] # ∂N/∂x
            dNdx[i, 2] = invJ[2, 1]*dNdxi[i, 1] + invJ[2, 2]*dNdxi[i, 2] # ∂N/∂y
        end

        # 5.4 Формирование матрицы B для элемента
        B = zeros(3, 8)
        for i in 1:4
            Bi = B_matrix(dNdx, i) 
            B[:, 2*i-1 : 2*i] .= Bi
        end

        # 5.5 Вклад текущей точки интегрирования в матрицу жесткости элемента
        BTDB = thickness * B' * D * B * detJ * weight
        element_stiffness[e] += BTDB
    end
end

# ***** 6. СОБОРКА ГЛОБАЛЬНОЙ МАТРИЦЫ ЖЕСТКОСТИ *****
K = assemble_global_stiffness(element_stiffness, connectivity)

# ***** 7. ЗАДАНИЕ ГРАНИЧНЫХ УСЛОВИЙ *****
fixed_dofs = Int[]


# ***** 8. ПРИЛОЖЕНИЕ НАГРУЗКИ (ДАВЛЕНИЯ) *****
# Инициализируем вектор внешних сил (по 2 степени свободы на узел)
F = zeros(2 * num_nodes)

# Расчет длины элемента по оси y
element_length_y = Ly / ny

# Сила на внутренний узел: давление * высота * толщина
node_force = pressure * element_length_y * thickness

# На крайние узлы (верхний и нижний) приходится половина силы
# так как они "обслуживают" только половину элемента
end_node_force = node_force / 2

# 1. Определяем закрепленные узлы (левая грань x=0)
fixed_dofs = Int[]
left_edge_nodes = []

for i in 1:num_nodes
    if abs(nodes[i, 1] - 0.0) < 1e-10
        push!(left_edge_nodes, i)
        push!(fixed_dofs, 2*i - 1)  # x-перемещение
        push!(fixed_dofs, 2*i)      # y-перемещение
    end
end

println("Закрепленные узлы на левой грани: $left_edge_nodes")

# 2. Приложение нагрузки (правая грань x=100)
F = zeros(2 * num_nodes)
right_edge_nodes = []

for i in 1:num_nodes
    if abs(nodes[i, 1] - 100.0) < 1e-10
        push!(right_edge_nodes, i)
        dof_x = 2*i - 1
        
        if abs(nodes[i, 2] - 0.0) < 1e-10 || abs(nodes[i, 2] - 50.0) < 1e-10
            F[dof_x] = end_node_force
        else
            F[dof_x] = node_force
        end
    end
end

println("Нагруженные узлы на правой грани: $right_edge_nodes")

# ********************** ФУНКЦИЯ ДЛЯ ВЫЧИСЛЕНИЯ ПОГРЕШНОСТИ В НОРМЕ L2
	function calculate_stress_error_L2(numerical_stresses, pressure)
    error_norm = 0.0
    exact_solution_norm = 0.0

    # Точное решение для задачи растяжения
    exact_sigma_xx = pressure
    exact_sigma_yy = 0.0
    exact_tau_xy = 0.0

    # Цикл по всем элементам
    for e in axes(numerical_stresses, 1)
        # Численное решение (напряжения в центре элемента)
        num_sigma_xx, num_sigma_yy, num_tau_xy = numerical_stresses[e]

        # Разности между численным и точным решениями
        diff_xx = num_sigma_xx - exact_sigma_xx
        diff_yy = num_sigma_yy - exact_sigma_yy
        diff_xy = num_tau_xy - exact_tau_xy

        # Вклад текущего элемента в норму ошибки
        error_norm += diff_xx^2 + diff_yy^2 + diff_xy^2

        # Вклад текущего элемента в норму точного решения
        exact_solution_norm += exact_sigma_xx^2 + exact_sigma_yy^2 + exact_tau_xy^2
    end

    # Относительная погрешность в норме L2
    if exact_solution_norm < 1e-12
        return (error_norm < 1e-12) ? 0.0 : 1.0
    else
        relative_error = sqrt(error_norm / exact_solution_norm)
        return relative_error
    end
end
	
# ***** 9. УЧЕТ ГРАНИЧНЫХ УСЛОВИЙ В СИСТЕМЕ *****
# Число степеней свободы всего
num_dofs = 2 * num_nodes

# Создаем карту DOF
function build_dof_map(num_dofs, fixed_dofs)
    dof_map = fill(-1, num_dofs)
    counter = 1

    for i in 1:num_dofs
        if !(i in fixed_dofs)
            dof_map[i] = counter
            counter += 1
        end
    end

    return dof_map
end

dof_map = build_dof_map(num_dofs, fixed_dofs)

# Инициализируем сокращенные матрицу жесткости и вектор сил
free_dofs = count(x -> x != -1, dof_map)  # число незакреплённых DOF
K_free = zeros(free_dofs, free_dofs)
F_free = zeros(free_dofs)

# Переносим данные из глобальной системы
for i in 1:num_dofs
    mi = dof_map[i]
    if mi == -1
        continue
    end

    F_free[mi] = F[i]

    for j in 1:num_dofs
        mj = dof_map[j]
        if mj == -1
            continue
        end

        K_free[mi, mj] = K[i, j]
    end
end

# ***** 10. РЕШЕНИЕ СИСТЕМЫ ЛИНЕЙНЫХ УРАВНЕНИЙ *****
# Решаем систему уравнений: K_free * U_free = F_free
U_free = K_free \ F_free

# ***** 11. ВОССТАНОВЛЕНИЕ ПОЛНОГО ВЕКТОРА ПЕРЕМЕЩЕНИЙ *****
# Полный вектор перемещений (все степени свободы, включая закреплённые)
num_dofs = 2 * num_nodes
U = zeros(num_dofs)

# Заполняем только свободные DOF из решения
for i in 1:num_dofs
    mi = dof_map[i]
    if mi != -1
        U[i] = U_free[mi]
    end
end


	# ***** 12. ВЫВОД РЕЗУЛЬТАТОВ *****

# 12.1 Вывод перемещений узлов
println("Node displacements (x, y):")
for i in 0:num_nodes-1
    x = nodes[i+1, 1]  # Julia: индексация с 1
    y = nodes[i+1, 2]
    u_x = U[2*i+1]
    u_y = U[2*i+2]
    println("Node $i at ($x, $y): (u_x = $u_x, u_y = $u_y)")
end

# 12.2 Расчет напряжений в элементах
element_stresses = calculate_stresses(element_stiffness, connectivity, U, nodes)


# 12.3 Вывод напряжений по элементам
println("\nElement stresses (sigma_xx, sigma_yy, tau_xy):")
for e in 1:length(element_stresses)
    # Вычисляем координаты центра элемента
    center_x = 0.0
    center_y = 0.0
    for i in 1:4
        node_id = connectivity[e, i]
        center_x += nodes[node_id, 1]
        center_y += nodes[node_id, 2]
    end
    center_x /= 4
    center_y /= 4

    # Получаем напряжения для элемента
    stress = element_stresses[e]

    # Выводим напряжения
    σxx = stress[1]
    σyy = stress[2]
    τxy = stress[3]

    println("Element $(e-1) at ($center_x, $center_y): [σxx = $σxx, σyy = $σyy, τxy = $τxy] МПа")
end

	# Вычисляем относительную погрешность
rel_error = calculate_stress_error_L2(element_stresses, pressure)

println("Relative L2 error: $(rel_error * 100)%")

	using LinearAlgebra

	# ошибки напряжений при разных точках интегрирования
function calculate_L2_error_gauss(numerical_stresses, connectivity, nodes, exact_solution; integration_order=2)
    # 1. Подготовка точек и весов Гаусса
    if integration_order == 1
        gauss_points = [0.0]
        gauss_weights = [2.0]
    elseif integration_order == 2
        gauss_points = [-1 / sqrt(3), 1 / sqrt(3)]
        gauss_weights = [1.0, 1.0]
    elseif integration_order == 3
        gauss_points = [-sqrt(3/5), 0.0, sqrt(3/5)]
        gauss_weights = [5/9, 8/9, 5/9]
    else
        error("Порядок интегрирования должен быть 1, 2 или 3")
    end

    # Общая норма ошибки и точного решения
    error_norm = 0.0
    exact_norm = 0.0

    # 2. Цикл по всем элементам
    for e in axes(connectivity, 1)
        node_coords = zeros(4, 2)
        for i in 1:4
            node_id = connectivity[e, i]
            node_coords[i, 1] = nodes[node_id, 1]  # x-координата узла
            node_coords[i, 2] = nodes[node_id, 2]  # y-координата узла
        end

        num_sol = numerical_stresses[e]  # численные напряжения для элемента

        # 3. Интегрирование по Гауссу (цикл по всем точкам)
        for (gp_xi, xi) in enumerate(gauss_points), (gp_eta, eta) in enumerate(gauss_points)
            weight = gauss_weights[gp_xi] * gauss_weights[gp_eta]

            # 3.2 Якобиан преобразования
            dNdxi = shape_function_derivatives(xi, eta)
            J = jacobian(dNdxi, node_coords)
            detJ = det(J)

            # 3.3 Вычисляем глобальные координаты точки интегрирования
            N = shape_functions(xi, eta)
            x = sum(N[i] * node_coords[i, 1] for i in 1:4)
            y = sum(N[i] * node_coords[i, 2] for i in 1:4)

            # 3.4 Точное решение в точке (x, y)
            exact_sol = exact_solution(x, y)

            # 3.5 Вклад в нормы
            diff = num_sol .- exact_sol
            error_norm += dot(diff, diff) * weight * abs(detJ)
            exact_norm += dot(exact_sol, exact_sol) * weight * abs(detJ)
        end
    end

    # 4. Относительная погрешность с защитой от деления на ноль
    eps = 1e-12
    if exact_norm < eps
        return error_norm < eps ? 0.0 : 1.0
    else
        return sqrt(error_norm / exact_norm)
    end
end
	# Аналитическое решение σxx = давление, σyy = 0, τxy = 0
exact_solution(x, y) = [pressure, 0.0, 0.0]
	rel_error = calculate_L2_error_gauss(element_stresses, connectivity, nodes, exact_solution, integration_order=1)

println("Относительная погрешность L₂: $(rel_error * 100)%")

	# ****** Функция, которая считает отн. погрешность только по компоненте σxx
	using LinearAlgebra

function calculate_sigma_xx_error_L2(numerical_stresses, connectivity, nodes, exact_sigma_xx; integration_order=2)
    # 1. Подготовка точек и весов Гаусса
    if integration_order == 1
        gauss_points = [0.0]
        gauss_weights = [2.0]
    elseif integration_order == 2
        gauss_points = [-1 / sqrt(3), 1 / sqrt(3)]
        gauss_weights = [1.0, 1.0]
    elseif integration_order == 3
        gauss_points = [-sqrt(3/5), 0.0, sqrt(3/5)]
        gauss_weights = [5/9, 8/9, 5/9]
    else
        error("Порядок интегрирования должен быть 1, 2 или 3")
    end

    error_sq = 0.0   # ∫(σxx_num - σxx_exact)^2 dΩ
    exact_sq = 0.0  # ∫(σxx_exact)^2 dΩ

    # 2. Цикл по всем элементам
    for e in axes(connectivity, 1)
        # 2.1. Координаты узлов элемента
        node_coords = zeros(4, 2)
        for i in 1:4
            node_id = connectivity[e, i]
            node_coords[i, 1] = nodes[node_id, 1]  # x-координата
            node_coords[i, 2] = nodes[node_id, 2]  # y-координата
        end

        # 2.2. Численное значение σxx в центре элемента (или среднее?)
        num_sigma_xx = numerical_stresses[e][1]

        # 3. Интегрирование по Гауссу
        for (gp_xi, xi) in enumerate(gauss_points), (gp_eta, eta) in enumerate(gauss_weights)
            weight = gauss_weights[gp_xi] * gauss_weights[gp_eta]

            # 3.2 Якобиан преобразования
            dNdxi = shape_function_derivatives(xi, eta)
            J = jacobian(dNdxi, node_coords)
            detJ = abs(det(J))

            # 3.3 Вычисляем глобальные координаты точки интегрирования
            N = shape_functions(xi, eta)
            x = sum(N[i] * node_coords[i, 1] for i in 1:4)
            y = sum(N[i] * node_coords[i, 2] for i in 1:4)

            # 3.4 Точное решение для σxx в этой точке
            exact_sigma = exact_sigma_xx(x, y)

            # 3.5 Вклад в интегралы
            diff = num_sigma_xx - exact_sigma
            error_sq += diff^2 * weight * detJ
            exact_sq += exact_sigma^2 * weight * detJ
        end
    end

    # 4. Относительная погрешность с защитой от деления на ноль
    eps = 1e-12
    if exact_sq < eps
        return error_sq < eps ? 0.0 : 1.0
    else
        return sqrt(error_sq / exact_sq)
    end
end
	exact_sigma_xx(x, y) = pressure  # постоянное напряжение вдоль оси x
error_xx = calculate_sigma_xx_error_L2(element_stresses, connectivity, nodes, exact_sigma_xx, integration_order=1)
println("Относительная погрешность σxx: $(error_xx * 100)%")

	using Plots
gr()

using Plots
gr()

	function calculate_von_mises_stress(element_stresses)
    von_mises = zeros(size(element_stresses, 1))
    for e in 1:size(element_stresses, 1)
        σxx, σyy, τxy = element_stresses[e]
        von_mises[e] = sqrt(σxx^2 - σxx*σyy + σyy^2 + 3*τxy^2)
    end
    return von_mises
end

function calculate_element_energy(element_stiffness, element_stresses, connectivity, nodes, U, thickness)
    num_elements = size(connectivity, 1)
    element_energy = zeros(num_elements)

    # Точки Гаусса и веса
    gauss_points = [-1 / sqrt(3), 1 / sqrt(3)]
    gauss_weights = [1.0, 1.0]

    for e in 1:num_elements
        node_ids = connectivity[e, :]

        # Получаем координаты узлов текущего элемента
        coords = zeros(4, 2)
        for i in 1:4
            coords[i, :] .= nodes[node_ids[i], :]
        end

        energy = 0.0

        for gp_xi in eachindex(gauss_points), gp_eta in eachindex(gauss_weights)
            xi = gauss_points[gp_xi]
            eta = gauss_points[gp_eta]
            weight = gauss_weights[gp_xi] * gauss_weights[gp_eta]

            # Производные функций формы
            dNdxi = shape_function_derivatives(xi, eta)

            # Якобиан
            J = jacobian(dNdxi, coords)
            detJ = det(J)
            invJ = inv(J)

            # Пересчет производных в глобальные координаты
            dNdx = zeros(4, 2)
            for i in 1:4
                dNdx[i, 1] = invJ[1, 1]*dNdxi[i, 1] + invJ[1, 2]*dNdxi[i, 2]
                dNdx[i, 2] = invJ[2, 1]*dNdxi[i, 1] + invJ[2, 2]*dNdxi[i, 2]
            end

            # Формирование матрицы B
            B = zeros(3, 8)
            for i in 1:4
                Bi = B_matrix(dNdx, i)
                B[:, 2*i - 1 : 2*i] .= Bi
            end

            # Перемещения узлов элемента
            u = zeros(8)
            for i in 1:4
                node = node_ids[i]
                u[2*i - 1] = U[2*node - 1]
                u[2*i] = U[2*node]
            end

            # Деформации: ε = B * u
            strains = B * u

            # Напряжения: σ = D * ε
            D = get_D_matrix()
            stresses = D * strains

            # Плотность энергии: U = ½ σᵀ ε
            density = 0.5 * dot(stresses, strains)

            # Вклад в энергию элемента
            energy += density * weight * abs(detJ) * thickness
        end

        element_energy[e] = energy
    end

    return element_energy
end

function plot_combined_field(element_stresses, connectivity, nodes, U; 
                              filename="output.png", 
                              plot_type=1,       # 1–3: σxx–τxy; 4–6: εxx–γxy; 7–8: u_x–u_y; 9: |u|
                              deform_scale=1.0,  # Масштаб деформации
                              deformed=false)    # Рисовать ли на деформированной сетке
    num_nodes = size(nodes, 1)
    num_elements = size(connectivity, 1)

    # Если требуется деформировать сетку — применяем перемещения
    if deformed
        new_nodes = copy(nodes)
        for i in 1:num_nodes
            u_x = U[2*i - 1]
            u_y = U[2*i]
            new_nodes[i, 1] += deform_scale * u_x
            new_nodes[i, 2] += deform_scale * u_y
        end
    else
        new_nodes = nodes
    end

    plt = plot(
               xlabel="x (мм)", ylabel="y (мм)",
               aspect_ratio=:equal, legend=false, dpi=300)

    cmap = palette(:viridis)

    # --- Отрисовка напряжений ---
    if plot_type in [1, 2, 3]
        stress_values = [element_stresses[e][plot_type] for e in 1:num_elements]
        vmin, vmax = extrema(stress_values)
        centers_x = zeros(num_elements)
        centers_y = zeros(num_elements)

        for e in 1:num_elements
            node_ids = connectivity[e, :]
            coords = [[new_nodes[i, 1], new_nodes[i, 2]] for i in node_ids]
            center_x = mean([c[1] for c in coords])
            center_y = mean([c[2] for c in coords])
            centers_x[e] = center_x
            centers_y[e] = center_y

            x = [coords[1][1], coords[2][1], coords[3][1], coords[4][1], coords[1][1]]
            y = [coords[1][2], coords[2][2], coords[3][2], coords[4][2], coords[1][2]]

            color_idx = clamp(Int(round((stress_values[e] - vmin) / (vmax - vmin + 1e-12) * (length(cmap) - 1))) + 1, 1, length(cmap))
            color = cmap[color_idx]

            plot!(x, y, fill=(true, 0.8, color), linecolor=:black, linewidth=0.5, label="")
        end

        scatter!(centers_x, centers_y, zcolor=stress_values, marker=(:square, 0),
                 clims=(vmin, vmax), color=:viridis,
                 colorbar_title=["    ", "    ", "    "][plot_type],
                 colorbar=true)

    # --- Отрисовка деформаций ---
    elseif plot_type in [4, 5, 6]
        component_index = plot_type - 3  # 4→1, 5→2, 6→3 для εxx, εyy, γxy
        D = get_D_matrix()
        invD = inv(D)  # Обратная матрица упругости

        strain_values = zeros(num_elements)
        for e in 1:num_elements
            σ = element_stresses[e]
            ε = invD * σ  # Вычисляем деформации из напряжений
            strain_values[e] = ε[component_index]
        end

        vmin, vmax = extrema(strain_values)

        for e in 1:num_elements
            node_ids = connectivity[e, :]
            coords = [[new_nodes[i, 1], new_nodes[i, 2]] for i in node_ids]
            x = [coords[1][1], coords[2][1], coords[3][1], coords[4][1], coords[1][1]]
            y = [coords[1][2], coords[2][2], coords[3][2], coords[4][2], coords[1][2]]
            s_avg = mean(strain_values[e])

            color_idx = clamp(Int(round((s_avg - vmin) / (vmax - vmin + 1e-12) * (length(cmap) - 1))) + 1, 1, length(cmap))
            color = cmap[color_idx]

            plot!(x, y, fill=(true, 0.8, color), linecolor=:black, linewidth=0.5, label="")
        end

        scatter!(zeros(0), zeros(0), zcolor=[vmin, vmax], c=:viridis, clims=(vmin, vmax),
                 colorbar_title=["    ", "    ", "    "][component_index],
                 colorbar=true)


    # --- Отрисовка компонент перемещений (u_x, u_y) ---
    elseif plot_type in [7, 8]
        component_index = plot_type == 7 ? 1 : 2  # 1: u_x, 2: u_y
        displacements = zeros(num_nodes)
        for i in 1:num_nodes
            displacements[i] = U[2*i - (component_index == 1 ? 1 : 0)]
        end
        vmin, vmax = extrema(displacements)

        for e in 1:num_elements
            node_ids = connectivity[e, :]
            coords = [[new_nodes[i, 1], new_nodes[i, 2]] for i in node_ids]
            x = [coords[1][1], coords[2][1], coords[3][1], coords[4][1], coords[1][1]]
            y = [coords[1][2], coords[2][2], coords[3][2], coords[4][2], coords[1][2]]
            disp_avg = mean(displacements[node_ids])

            color_idx = clamp(Int(round((disp_avg - vmin) / (vmax - vmin + 1e-12) * (length(cmap) - 1))) + 1, 1, length(cmap))
            color = cmap[color_idx]

            plot!(x, y, fill=(true, 0.8, color), linecolor=:black, linewidth=0.5, label="")
        end

        scatter!(zeros(0), zeros(0), zcolor=[vmin, vmax], c=:viridis, clims=(vmin, vmax),
                 colorbar_title=plot_type == 7 ? "    " : "    ",
                 colorbar=true)

    # --- Отрисовка модуля перемещений ---
    elseif plot_type == 9
        displacements = zeros(num_nodes)
        for i in 1:num_nodes
            u_x = U[2*i - 1]
            u_y = U[2*i]
            displacements[i] = sqrt(u_x^2 + u_y^2)
        end
        vmin, vmax = extrema(displacements)

        for e in 1:num_elements
            node_ids = connectivity[e, :]
            coords = [[new_nodes[i, 1], new_nodes[i, 2]] for i in node_ids]
            x = [coords[1][1], coords[2][1], coords[3][1], coords[4][1], coords[1][1]]
            y = [coords[1][2], coords[2][2], coords[3][2], coords[4][2], coords[1][2]]
            disp_avg = mean(displacements[node_ids])

            color_idx = clamp(Int(round((disp_avg - vmin) / (vmax - vmin + 1e-12) * (length(cmap) - 1))) + 1, 1, length(cmap))
            color = cmap[color_idx]

            plot!(x, y, fill=(true, 0.8, color), linecolor=:black, linewidth=0.5, label="")
        end

        scatter!(zeros(0), zeros(0), zcolor=[vmin, vmax], c=:viridis, clims=(vmin, vmax),
                 colorbar_title="   ", colorbar=true)

		elseif plot_type == 10
    element_energy = calculate_element_energy(element_stiffness, element_stresses, connectivity, nodes, U, thickness)
    vmin, vmax = extrema(element_energy)

    for e in 1:num_elements
        node_ids = connectivity[e, :]
        coords = [[new_nodes[i, 1], new_nodes[i, 2]] for i in node_ids]
        x = [coords[1][1], coords[2][1], coords[3][1], coords[4][1], coords[1][1]]
        y = [coords[1][2], coords[2][2], coords[3][2], coords[4][2], coords[1][2]]

        color_idx = clamp(Int(round((element_energy[e] - vmin) / (vmax - vmin + 1e-12) * (length(cmap) - 1))) + 1, 1, length(cmap))
        color = cmap[color_idx]

        plot!(x, y, fill=(true, 0.8, color), linecolor=:black, linewidth=0.5, label="")
    end

    scatter!(zeros(0), zeros(0), zcolor=[vmin, vmax], c=:viridis, clims=(vmin, vmax),
             colorbar_title="    ", colorbar=true)

		# --- Отрисовка напряжений Мизеса ---
elseif plot_type == 11
    von_mises = calculate_von_mises_stress(element_stresses)
    vmin, vmax = extrema(von_mises)

    for e in 1:num_elements
        node_ids = connectivity[e, :]
        coords = [[new_nodes[i, 1], new_nodes[i, 2]] for i in node_ids]
        x = [coords[1][1], coords[2][1], coords[3][1], coords[4][1], coords[1][1]]
        y = [coords[1][2], coords[2][2], coords[3][2], coords[4][2], coords[1][2]]
        vm_avg = mean(von_mises[e])

        color_idx = clamp(Int(round((vm_avg - vmin) / (vmax - vmin + 1e-12) * (length(cmap) - 1))) + 1, 1, length(cmap))
        color = cmap[color_idx]

        plot!(x, y, fill=(true, 0.8, color), linecolor=:black, linewidth=0.5, label="")
    end

    scatter!(zeros(0), zeros(0), zcolor=[vmin, vmax], c=:viridis, clims=(vmin, vmax),
             colorbar_title="    ", colorbar=true)
    else
        error("Неверный тип графика. Используйте:\n" *
              "plot_type=1 (σxx), 2 (σyy), 3 (τxy)\n" *
              "plot_type=4 (εxx), 5 (εyy), 6 (γxy)\n" *
              "plot_type=7 (u_x), 8 (u_y), 9 (|u|)")
    end

    display(plt)
    savefig(plt, filename)
end

	# Отрисовка σxx
plot_combined_field(element_stresses, connectivity, nodes, U, 
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\sigma_xx_combined.png", 
                    plot_type=1)

# Отрисовка σyy
plot_combined_field(element_stresses, connectivity, nodes, U, 
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\sigma_yy_combined.png", 
                    plot_type=2)

# Отрисовка τxy
plot_combined_field(element_stresses, connectivity, nodes, U, 
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\tau_xy_combined.png", 
                    plot_type=3)

# Отрисовка модуля перемещений
plot_combined_field(element_stresses, connectivity, nodes, U, 
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\displacement_combined.png", 
                    plot_type=9)

	# Отрисовка u_x
plot_combined_field(element_stresses, connectivity, nodes, U,
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\u_x.png",
                    plot_type=7)

# Отрисовка u_y
plot_combined_field(element_stresses, connectivity, nodes, U,
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\u_y.png",
                    plot_type=8)
# Деформации
plot_combined_field(element_stresses, connectivity, nodes, U,
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\strain_exx.png",
                    plot_type=4)

plot_combined_field(element_stresses, connectivity, nodes, U,
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\strain_eyy.png",
                    plot_type=5)

plot_combined_field(element_stresses, connectivity, nodes, U,
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\strain_gxy.png",
                    plot_type=6)

	# Вычисляем поэлементную энергию
element_energy = calculate_element_energy(element_stiffness, element_stresses, connectivity, nodes, U, thickness)

# Отрисовка энергии
plot_combined_field(element_stresses, connectivity, nodes, U;
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\energy_distribution.png",
                    plot_type=10)
	# Отрисовка напряжения Мизеса
plot_combined_field(element_stresses, connectivity, nodes, U,
                    filename="C:\\Users\\Admin\\Music\\RESULTS\\von_mises.png",
                    plot_type=11)

# ******** ПОДСЧЕТ ПОЛНОЙ ЭНЕРГИИ ********
total_energy = sum(calculate_element_energy(element_stiffness, element_stresses, connectivity, nodes, U, thickness))
println("Общая энергия упругой деформации: $total_energy МПа·мм³")
inv_dx = nx / Lx
println("1/шаг по оси x = $inv_dx 1/мм")

	# ******** НАПРЯЖЕНИЕ ПО ХХ отдельно ******
function calculate_mean_sigma_xx(element_stresses)
    num_elements = length(element_stresses)
    total = 0.0
    for e in 1:num_elements
        total += element_stresses[e][1]
    end
    return total / num_elements
end

	# Расчёт среднего σxx
mean_sigma_xx = calculate_mean_sigma_xx(element_stresses)

# Вывод результата
println("Среднее значение σxx: $mean_sigma_xx МПа")

end


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Plots = "~1.40.19"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "f2ebe675cefcceb97e5403bea52949af13b1648d"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "6c72198e6a101cccdd4c9731d3985e904ba26037"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7bb1361afdb33c7f2b085aa49ea8fe1b0fb14e58"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.1+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "83dc665d0312b41367b7263e8a4d172eac1897f4"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.4"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3a948313e7a41eb1db7a1e733e6335f17b4ab3c4"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "7.1.1+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "1828eb7275491981fa5f1752a5e126e8f26f8741"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.17"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "27299071cc29e409488ada41ec7643e0ab19091f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.17+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "35fbd0cefb04a516104b8e183ce0df11b70a3f1a"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e95866623950267c1e4878846f848d94810de475"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "706dfd3c0dd56ca090e86884db6eda70fa7dd4af"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d3c8af829abaeba27181db4acb485b18d15d89c6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2ae7d4ddec2e13ad3bddf5c0796f7547cf682391"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.2+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "275a9a6d85dc86c24d03d1837a0010226a96f540"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.3+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "0c5a5b7e440c008fe31416a3ac9e0d2057c81106"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.19"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "eb38d376097f47316fe089fc62cb7c6d85383a52"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "e1d5e16d0f65762396f9ca4644a5f4ddab8d452b"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2c962245732371acd51700dbb268af311bddd719"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.6"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "6258d453843c466d84c17a58732dda5deeb8d3af"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.24.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "af305cc62419f9bd61b6644d19170a4d258c7967"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.7.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "9caba99d38404b285db8801d5c45ef4f4f425a6d"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.1+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "c5bf2dad6a03dfef57ea0a170a1fe493601603f2"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.5+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4bba74fa59ab0755167ad24f98800fe5d727175b"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.12.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "fbf139bce07a534df0e699dbb5f5cc9346f95cc1"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.9.2+0"
"""

# ╔═╡ Cell order:
# ╠═9bc01b20-8bc2-11f0-32ac-a587080cee1e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
