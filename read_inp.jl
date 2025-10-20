begin

	####################### ЧИТКА ИЗ ФАЙЛА ##############################
	using DelimitedFiles
	using LinearAlgebra
	using Statistics 
	
	function read_inp_file(filename)
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
	

	Lx, Ly, nx, ny, nodes, connectivity = read_inp_file("C:\\Users\\Admin\\res2.inp")
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
