from sympy import sympify,expand,factor,Mul,solve,Symbol
print ("您正在使用Marshall矩阵计算机.You're using Marshall matrix computer\n按行即按1~n行顺序从左往右输入，数字间只空一格，例：\n|1 2|\n|3 4| 输入：'1 2 3 4'")
line = '--'*60
print(line)
import fractions

def clean(m):
    l = m.split(' ')
    for i in range(len(l)):
        if '**' in l[i]:
            p = l[i].split('**')
            l[i] = float(p[0]) ** float(p[1])
            l[i] = str(l[i])
        if '*' in l[i]:
            p = l[i].split('*')
            l[i] = float(p[0]) * float(p[1])
    return list(map(float, l))

def stairs_det(n):
    if n%1!=0 or n<0:
        return False
    elif n==0:
        return 1
    return n*stairs_det(n-1)

def stairs_m(n):
    if n%1!=0 or n<0:
        return False
    elif n==0:
        return 1
    return n**2*stairs_m(n-1)

def determinant_laplace(det):
    l = len(det)**(1/2)

    def cal_2(det):
        return det[0]*det[3]-det[2]*det[1]

    if l%1!=0:
        print('\n非方阵')
    elif l == 1:
        print('1阶\n行列式值为：',det[0])
    else:
        l = int(l)
        print(f'{l}阶')

        v = {}
        for i in range(0,l**2,l):
            #将每行单列一个列表，放入字典
            v_name = f'row{i//l}'
            v[v_name] = det[i:i+l]

        for row,val in v.items():
            print(f'{row} : {list(map(str,val))}')

        if l>9:
            print(f'\n{l}阶行列式需要计算{stairs_det(l)}次，请您稍安勿躁...')

        def expansion(det,l):
            a_n = []
            first_row = {}
            row_name = f'row_1st{l}'
            first_row[row_name] = det[0:l]
            son_dic = {}
            if l<3:
                return cal_2(det)

            for i in range(l):
                son_name = f'son{l,i}'
                son_dic[son_name] = det[l:l ** 2]#取l个去除原式第一行的向量
            for i in range(l):
                for j in range(l-1):
                    del son_dic[f'son{l,i}'][i+j*l-j]#可以得到任意阶数行列式以第一行展开的子式

            if l>2:
                for i in range(l):
                    a_n.append((-1)**i*(first_row[f'row_1st{l}'][i])*(expansion(son_dic[f'son{l,i}'],l-1)))#递归到二阶
            return sum(a_n)
        print('\n行列式值为:',expansion(det,l))

def determinant_gauss(det):
    l = len(det) ** (1 / 2)
    if l%1!=0:
        print('\n非方阵')
    elif l == 1:
        print('1阶\n行列式值为：',det[0])
    else:
        l = int(l)
        v = {}
        for i in range(0, l ** 2, l):
            # 将每行单列一个列表，放入字典
            v_name = f'row{i // l}'
            v[v_name] = det[i:i + l]

        for row,val in v.items():
            print(f'{row} : {list(map(str,val))}')

        diagonal = []
        quotient = {}
        for i in range(l):
            for j in range(i+1,l):
                for k in range(i,l):
                    if v[f'row{i}'][i] == 0:#如果(i,i)位置为0，则将下一行换上来
                        if k==l-1:
                            print('\n行列式值为：0')#换到l-1还是0，说明全是0
                            return 0
                        v[f'row{i}'],v[f'row{k+1}'] = v[f'row{k+1}'],v[f'row{i}']
                    else:
                        break
                quotient[f'q{i, j}'] = v[f'row{j}'][i]/v[f'row{i}'][i]#把j行的第一个元素与a(i,i)的比求出来

                for k in range(l):
                    v[f'row{j}'][k] -= quotient[f'q{i, j}']*v[f'row{i}'][k]#每行逐个相减
            diagonal.append(v[f'row{i}'][i])

        def multi(list, n):  # 从列表第n个乘到第0个
            if n < 0:
                return 1
            return list[n] * multi(list, n - 1)
        print(f'\n行列式值为：{multi(diagonal,l-1)}')

def inverse(det):
    l = len(det) ** (1 / 2)
    if l%1 != 0 :
        print('非方阵，无法计算')
    elif l ==1:
        if det[0] == 0:
            print('行列式为零，无逆矩阵')
        else:
            print(f'逆矩阵为：{det[0]**(-1)}')
    else:
        def cal_2(det):
            if len(det)==1:
                return det[0]
            else:
                return det[0] * det[3] - det[2] * det[1]
        l = int(l)
        print(f'{l}阶')

        v = {}
        for i in range(0, l ** 2, l):
            v_name = f'row{i // l}'
            v[v_name] = det[i:i + l]

        for row,ele in v.items():
            print(row, '|', " | ".join([f"{str(val):>8}" for val in ele]) + " |")

        def expansion(det, l):
            a_n = []
            first_row = {}
            row_name = f'row_1st{l}'
            first_row[row_name] = det[0:l]
            son_dic = {}
            if l < 3:
                return cal_2(det)

            for i in range(l):
                son_name = f'son{l, i}'
                son_dic[son_name] = det[l:l ** 2]  # 取l个去除原式第一行的向量
            for i in range(l):
                for j in range(l - 1):
                    del son_dic[f'son{l, i}'][i + j * l - j]  # 可以得到任意阶数行列式以第一行展开的子式

            if l > 2:
                for i in range(l):
                    a_n.append((-1) ** i * (first_row[f'row_1st{l}'][i]) * (expansion(son_dic[f'son{l, i}'], l - 1)))
            return sum(a_n)

        def determinant_gauss(det,l):

            v = {}
            for i in range(0, l ** 2, l):
                # 将每行单列一个列表，放入字典
                v_name = f'row{i // l}'
                v[v_name] = det[i:i + l]


            diagonal = []
            quotient = {}
            for i in range(l):
                for j in range(i + 1, l):
                    for k in range(i, l):
                        if v[f'row{i}'][i] == 0:  # 如果(i,i)位置为0，则将下一行换上来
                            if k == l - 1:
                                return 0  # 换到l-1还是0，说明全是0
                            v[f'row{i}'], v[f'row{k + 1}'] = v[f'row{k + 1}'], v[f'row{i}']
                        else:
                            break
                    quotient[f'q{i, j}'] = v[f'row{j}'][i] / v[f'row{i}'][i]  # 把j行的第一个元素与a(i,i)的比求出来

                    for k in range(l):
                        v[f'row{j}'][k] -= quotient[f'q{i, j}'] * v[f'row{i}'][k]  # 每行逐个相减
                diagonal.append(v[f'row{i}'][i])

            def multi(list, n):  # 从列表第n个乘到第0个
                if n < 0:
                    return 1
                return list[n] * multi(list, n - 1)
            return multi(diagonal,l-1)

        if l<9:
            det_val = expansion(det, l)
            print('\n行列式值为:', expansion(det, l))
        else:
            det_val = determinant_gauss(det, l)
            print('\n行列式值为:', determinant_gauss(det,l))
        if expansion(det, l)  == 0 or determinant_gauss(det,l) == 0 :
            print('无逆矩阵')
        else:#只判断前一个if是否为真
            son_dic = {}
            for i in range(l):
                for j in range(l):
                    son_name = f'son{i,j}'
                    son_dic[son_name] = det[0 :l**2]
                    del son_dic[f'son{i,j}'][i*l:(i+1)*l]
            for i in range(l**2):
                for j in range(l - 1):
                    del son_dic[f'son{i//l,i-i//l*l}'][i-(i//l*l)+j*l-j]#得到所有元素子式
            inverse_matrix ={}
            for i in range(l):
                in_name = f'inversed_row{i}'
                inverse_matrix[in_name] = []

            for i in range(l):
                for j in range(l):
                    if l<9:
                        inverse_matrix[f'inversed_row{j}'].append(((-1)**(i+j)*expansion(son_dic[f'son{i,j}'], l-1))/det_val)
            print('\n逆矩阵为：')
            for row,ele in inverse_matrix.items():
                print(row, '|', " | ".join([f"{str(val):>8}" for val in ele]) + " |")

def adjoint(det):
    l = len(det) ** (1 / 2)
    if l % 1 != 0:
        print('非方阵，无法计算')
    elif l == 1:
        if det[0] == 0:
            print('伴随矩阵为0或无定义')
        else:
            print('伴随矩阵为单位矩阵或无定义')
    else:
        def cal_2(det):
            if len(det) == 1:
                return det[0]
            else:
                return det[0] * det[3] - det[2] * det[1]

        l = int(l)
        print(f'{l}阶')

        v = {}
        for i in range(0, l ** 2, l):
            v_name = f'row{i // l}'
            v[v_name] = det[i:i + l]

        def expansion(det, l):
            a_n = []
            first_row = {}
            row_name = f'row_1st{l}'
            first_row[row_name] = det[0:l]
            son_dic = {}
            if l < 3:
                return cal_2(det)

            for i in range(l):
                son_name = f'son{l, i}'
                son_dic[son_name] = det[l:l ** 2]  # 取l个去除原式第一行的向量
            for i in range(l):
                for j in range(l - 1):
                    del son_dic[f'son{l, i}'][i + j * l - j]  # 可以得到任意阶数行列式以第一行展开的子式

            if l > 2:
                for i in range(l):
                    a_n.append(
                        (-1) ** i * (first_row[f'row_1st{l}'][i]) * (expansion(son_dic[f'son{l, i}'], l - 1)))
            return sum(a_n)
        def determinant_gauss(det,l):

            v = {}
            for i in range(0, l ** 2, l):
                # 将每行单列一个列表，放入字典
                v_name = f'row{i // l}'
                v[v_name] = det[i:i + l]

            diagonal = []
            quotient = {}
            for i in range(l):
                for j in range(i + 1, l):
                    for k in range(i, l):
                        if v[f'row{i}'][i] == 0:  # 如果(i,i)位置为0，则将下一行换上来
                            if k == l - 1:
                                return 0  # 换到l-1还是0，说明全是0
                            v[f'row{i}'], v[f'row{k + 1}'] = v[f'row{k + 1}'], v[f'row{i}']
                        else:
                            break
                    quotient[f'q{i, j}'] = v[f'row{j}'][i] / v[f'row{i}'][i]  # 把j行的第一个元素与a(i,i)的比求出来

                    for k in range(l):
                        v[f'row{j}'][k] -= quotient[f'q{i, j}'] * v[f'row{i}'][k]  # 每行逐个相减
                diagonal.append(v[f'row{i}'][i])

            def multi(list, n):  # 从列表第n个乘到第0个
                if n < 0:
                    return 1
                return list[n] * multi(list, n - 1)
            return multi(diagonal,l-1)

        son_dic = {}
        for i in range(l):
            for j in range(l):
                son_name = f'son{i, j}'
                son_dic[son_name] = det[0:l ** 2]
                del son_dic[f'son{i, j}'][i * l:(i + 1) * l]
        for i in range(l ** 2):
            for j in range(l - 1):
                del son_dic[f'son{i // l, i - i // l * l}'][i - (i // l * l) + j * l - j]
        adjoint_matrix = {}
        for i in range(l):
            ad_name = f'adjoint_row{i}'
            adjoint_matrix[ad_name] = []

        for i in range(l):
            for j in range(l):
                if l < 9:
                    adjoint_matrix[f'adjoint_row{i}'].append(
                        ((-1) ** (i + j) * expansion(son_dic[f'son{j, i}'], l - 1)))
                else:
                    adjoint_matrix[f'adjoint_row{i}'].append(((-1) ** (i + j) * determinant_gauss(son_dic[f'son{j, i}'],l-1)))

        print('\n伴随矩阵为：')
        for row, ele in adjoint_matrix.items():
            print(row,'|'," | ".join([f"{str(val):>8}" for val in ele]) + " |")

def matrix_multiply(det1,det2,t,power=False):
    mat_dic1 = {}
    mat_dic2 = {}
    for i in range(r1):
        mat_dic1[f'm1_row{i}'] = det1[i*c1:(i+1)*c1]
    mat_dic2r = {}
    for i in range(c2):
        mat_dic2[f'm2_col{i}'] = []
        for j in range(r2):
            mat_dic2[f'm2_col{i}'].append(det2[j*c2+i])
    for i in range(r2):
        mat_dic2r[f'm2_row{i}'] = det2[i*c2:(i+1)*c2]
    global Product_m
    Product_m = {}
    product = {}
    for i in range(r1):
        Product_m[f'product_m_row{i}'] = []
        for j in range(c2):
            product[f'product_row{j}'] = []
            for k in range(c1):
                product[f'product_row{j}'].append(mat_dic1[f'm1_row{i}'][k]*mat_dic2[f'm2_col{j}'][k])
            Product_m[f'product_m_row{i}'].append(sum(product[f'product_row{j}']))
    if not power:
        print('\n矩阵1')
        for row, ele in mat_dic1.items():
            print('       ',row, '|', " | ".join([f"{str(val):>8}" for val in ele]) + " |")
        print('\n矩阵2')
        for row, ele in mat_dic2r.items():
            print('       ',row, '|', " | ".join([f"{str(val):>8}" for val in ele]) + " |")
    if not power or t == times-1 :
        print('\n矩积为：')
        for row,lst in Product_m.items():
            print(row, '|', " | ".join([f"{str(val):>8}" for val in lst]) + " |")
    def det_product(Product_m,r1):
        det_p = []
        for i in range(r1):
            det_p += Product_m[f'product_m_row{i}']
        return det_p
    return det_product

def unknown_det(det,det_val):
    factor_bool = False
    l = len(det) ** (1 / 2)
    if l % 1 != 0:
        print('\n非方阵')
    elif l == 1:
        print('1阶\n行列式值为：', det[0])
    else:
        l = int(l)
        v = {}
        for i in range(0, l ** 2, l):
            # 将每行单列一个列表，放入字典
            v_name = f'row{i // l}'
            v[v_name] = det[i:i + l]
        for row, val in v.items():
            print(f'{row} : {list(map(str, val))}')
        element = {}


        def cal_def(v, l):
            for i in range(l):
                key = 'row' + str(i)
                for j in range(l):
                    e_name = f'a{i, j}'
                    element[e_name] = v[key][j]

            def cal_order():
                import itertools
                global factor_bool
                def order(l):
                    nums = list(range(l))
                    return list(itertools.permutations(nums))

                cal = {}
                for tup in order(l):
                    cal[tup] = []
                    for i in range(l):
                        cal[tup].append(element[f'a{i, tup[i]}'])

                consequence = 0

                def invs(l):
                    m = 0
                    n = len(l)
                    for i in range(n):
                        for j in range(i + 1, n):
                            if l[i] > l[j]:
                                m += 1
                    return m

                for name, val in cal.items():
                    val = sympify(val)
                    sign = (-1) ** invs(name)
                    consequence += Mul(sign,*val)
                expanded = expand(consequence)
                if det_val == 0:
                    print('行列式为：',expanded)
                    if expanded != factor(consequence):
                        factor_bool = True
                        print(f'因式分解：{factor(consequence)}')
                    else:
                        print('（无法因式分解）\n')
                else:
                    print(f'行列式：{expanded} = {det_val}')
                    if expanded-det_val != factor(consequence-det_val):
                        factor_bool = True
                        print(f'因式分解：{factor(consequence-det_val)}')
                    else:
                        print('（无法因式分解）\n')
                return expanded - det_val

            return cal_order
        '''
        def cal_one_unknown(func):
            exp = True
            try:
                expr = sympify(func)
                def f(x):
                    return float(eval(expr,{'x': x}))

            except Exception as e:
                exp = False
                print(f'Error: {e}，请输入正确的表达式！\n')
                def f(x):
                    return 1
            
            try:
                gap = float(input('请输入根之间的最小距离（大概；若仅有一根则输入该根粗略绝对值）：')) / 10
            except:
                print('（默认0.04）')
                gap = 0.04
            x1 = float(0)
            x2 = float(gap)
            x3 = x1
            x4 = x2
            t = 0
            r_t = 0
            t_1 = 0
            a = 1e-15
            a1 = 1e-10
            a2 = 1e-8
            try:
                r_num = int(input('请输入（大于或等于）根的数量：'))
            except:
                r_num = 3
                print('（默认3）')
            print("开始计算，请等待至显示'计算完毕'")
            if exp:
                for i in range(r_num + 1):  # 1根数
                    r_t += 1
                    t0 = 0
                    t = t - t

                    if f(x1) * f(x2) is not complex:

                        while abs(f((x1 + x2) / 2)) > a and abs(x1 - x2) > a:
                            mid = (x1 + x2) / 2
                            t0 += 1

                            if abs(x1 - x2) <= a and abs(f((x1 + x2) / 2)) > a2:
                                break
                            elif abs(f(mid)) > abs(f(x1)) and abs(f(mid)) > abs(f(x2)):
                                break
                            elif t0 >= 4000:  # 2精确度&逃脱循环&追找步数
                                break

                            elif abs(f(mid)) < abs(f(x2)) <= abs(f(x1)):
                                x1 = mid
                            elif abs(f(mid)) < abs(f(x1)) <= abs(f(x2)):
                                x2 = mid
                            elif abs(f(x1)) > abs(f(mid)) > abs(f(x2)):
                                x1 = mid
                                x2 = 2 * x2 - x1
                            else:
                                x2 = mid
                                x1 = 2 * x1 - x2

                        if abs(f((x1 + x2) / 2)) <= a1:
                            if f((x1 + x2) / 2) != 0:
                                print(f'根{r_t}=', (x1 + x2) / 2, '.', sep='')

                            elif f((x1 + x2) / 2) == 0:
                                print(f'根{r_t}=', (x1 + x2) / 2, '.', sep='')
                        else:
                            r_t = r_t - 1

                    x1 = x2
                    x2 = x1 + gap  # 3粗略度
                    while abs(x1 - x2) > a:
                        mid = (x1 + x2) / 2
                        t += 1

                        if t >= 8000:  # 同2
                            break
                        elif abs(f(mid)) > abs(f(x2)) >= abs(f(x1)):
                            x1 = mid
                        elif abs(f(mid)) > abs(f(x1)) >= abs(f(x2)):
                            x2 = mid
                        elif abs(f(x1)) < abs(f(mid)) < abs(f(x2)):
                            x1 = mid
                            x2 = 2 * x2 - x1
                        else:
                            x2 = mid
                            x1 = 2 * x1 - x2

                    else:
                        pass

                    x1 = x2
                    x2 = x1 + gap  # 同3

                for i in range(r_num):  # 1根数
                    r_t += 1
                    t0 = 0
                    t = t - t
                    t_1 += 1

                    if f(x3) * f(x4) is not complex:

                        while abs(f((x3 + x4) / 2)) > a and abs(x3 - x4) > a:
                            mid = (x3 + x4) / 2
                            t0 += 1
                            if abs(x3 - x4) <= a and abs(f((x3 + x4) / 2)) > a2:
                                break
                            elif abs(f(mid)) > abs(f(x3)) and abs(f(mid)) > abs(f(x4)):
                                break
                            elif t0 >= 4000:  # 2精确度&逃脱循环&追找步数
                                break
                            elif abs(f(mid)) < abs(f(x4)) <= abs(f(x3)):
                                x3 = mid
                            elif abs(f(mid)) < abs(f(x3)) <= abs(f(x4)):
                                x4 = mid
                            elif abs(f(x3)) > abs(f(mid)) > abs(f(x4)):
                                x3 = mid
                                x4 = 2 * x4 - x3
                            else:
                                x4 = mid
                                x3 = 2 * x3 - x4

                        if abs(f((x3 + x4) / 2)) <= a1 and t_1 >= 2:
                            print(f'根{r_t}=', (x3 + x4) / 2, '.', sep='')

                        elif f((x3 + x4) / 2) == 0 and t_1 >= 2:
                            print(f'根{r_t}=', (x3 + x4) / 2, '.', sep='')
                        else:
                            r_t = r_t - 1

                    x4 = x3
                    x3 = x4 - gap  # 3粗略度
                    while abs(x3 - x4) > a:
                        mid = (x3 + x4) / 2
                        t += 1
                        if abs(f(mid)) < abs(f(x3)) and abs(f(mid)) < abs(f(x4)):
                            break
                        elif t >= 8000:  # 同2
                            break
                        elif abs(f(mid)) > abs(f(x4)) >= abs(f(x3)):
                            x3 = mid
                        elif abs(f(mid)) > abs(f(x3)) >= abs(f(x4)):
                            x4 = mid
                        elif abs(f(x3)) < abs(f(mid)) < abs(f(x4)):
                            x3 = mid
                            x4 = 2 * x4 - x3
                        else:
                            x4 = mid
                            x3 = 2 * x3 - x4

                    else:
                        pass

                    x4 = x3
                    x3 = x4 - gap  # 同3
        '''
        f = cal_def(v, l)()
        skip = input('是否跳过计算（Y/n）：')
        if skip == 'n':
            x = Symbol('x')
            ans = solve(str(f),x)
            print(ans,'\n计算完毕，无更多解')

def fss(det):# Fundamental set of solutions
    l = len(det) ** (1 / 2)
    if l % 1 != 0:
        print('\n非矩阵或仅有零解，无基础解系')
    else:
        l = int(l)
        v = {}
        for i in range(0, l ** 2, l):
            # 将每行单列一个列表，放入字典
            v_name = f'row{i // l}'
            v[v_name] = det[i:i + l]
        print()
        for row, val in v.items():
            print(f'{row} : {list(map(str,val))}')
        print()

        pivot_rows = []
        quotient = {}
        fig = True
        seek_col = 0  #But failed
        pivot_col = []
        t = 0
        for i in range(l):
            pivot_bool = False
            i2 = i - seek_col
            for j in range(i2 + 1 if i<l-1 else i2, l):
                if not pivot_bool:
                    for k in range(i2, l):
                        if v[f'row{i2}'][i] == 0:  # 如果(i2,i)位置为0，则将下一行换上来 'Implementation of finding pivot'
                            if k == l - 1:
                                seek_col += 1
                                fig = False
                                break
                            v[f'row{i2}'], v[f'row{k + 1}'] = v[f'row{k + 1}'], v[f'row{i2}']
                        else:
                            t += 1
                            pivot_bool = True
                            pivot_rows.append(v[f'row{i2}'])
                            fig = True
                            pivot_col.append(seek_col+t-1)
                            break
                if fig and i != l - 1:
                    quotient[f'q{i2, j}'] = v[f'row{j}'][i] / v[f'row{i2}'][i]  # 把j行的第一个元素与a(i2,i)的比求出来

                    for k in range(l):
                        v[f'row{j}'][k] -= quotient[f'q{i2, j}'] * v[f'row{i2}'][k]  # 每行逐个相减
                else:
                    break
        R = l - seek_col
        for i in range(seek_col):
            pivot_rows.append([0]*l)
        for i,lst in enumerate(pivot_rows):
            print(f'echelon_row{i} : {list(map(str,lst))}')

        print()
        reduced_rows = {}
        for i in range(l):
            if i <R:
                reduced_rows[f'reduced_row{i}'] = [x/pivot_rows[i][pivot_col[i]] for x in pivot_rows[i]]
            else:
                reduced_rows[f'reduced_row{i}'] = ([0]*l)
            print(f'reduced_row{i} : {list(map(str,reduced_rows[f'reduced_row{i}']))}')

        print()
        standard_rows = reduced_rows.copy()
        for i in range(l):
            if i <R:
                for j in range(i+1,R):
                    standard_rows[f'reduced_row{i}'] = [a - b*standard_rows[f'reduced_row{i}'][j] for a,b in zip(standard_rows[f'reduced_row{i}'], standard_rows[f'reduced_row{j}'])]
            else:
                standard_rows[f'reduced_row{i}'] = ([0]*l)
            print(f'standard_row{i} : {list(map(str,standard_rows[f'reduced_row{i}']))}')

        print()
        if seek_col == 0:
            print(f'\nrank=n={R}，无基础解系')
        Xi = {}
        for i in range(seek_col):
            Xi[f'xi_{i}'] = []
            for j in range(l):
                Xi[f'xi_{i}'].append(-standard_rows[f'reduced_row{j}'][l-1-i])
            Xi[f'xi_{i}'][l-i-1] = 1
            print(f'xi_{i}:')
            print(" | \n".join([f"{str(val):>8}" for val in Xi[f'xi_{i}']]) + " |")

user_input = ''
if __name__ == '__main__':
    while True:
        USER_input = input('请选择（输入字母：a.行列式；b.伴随矩阵；c.逆矩阵；d.矩阵乘法；e.未知数矩阵；f.基础解系；g.特征向量）：').strip().lower()
        if USER_input =='q':
            break
        if USER_input =='a':
            while True:
                try:
                    det = input("请输入行列式：").strip()
                    if det == '':
                        break
                    det = list(map(str, det.split(' ')))
                    det = [fractions.Fraction(item) for item in det]
                    if len(det) <80:
                        determinant_laplace(det)
                    else:
                        determinant_gauss(det)
                    print('')
                    user_input = input("按回车重新运行，输入'q'退出，输入'b'回到首页：").strip().lower()
                    print(line)
                    if user_input=='q' or user_input=='b':
                        break
                except Exception as e:
                    print(f'Error: {e}:')
        if USER_input =='b':
            while True:
                try:
                    det = input("请输入矩阵：").strip()
                    if det =='':
                        break
                    det = list(map(str, det.split(' ')))
                    det = [fractions.Fraction(item) for item in det]
                    determinant_laplace(det)
                    adjoint(det)
                    print('')
                    user_input = input("按回车重新运行，输入'q'退出，输入'b'回到首页：").strip().lower()
                    print(line)
                    if user_input == 'q' or user_input == 'b':
                        break
                except Exception as e:
                    print(f'Error: {e}:')
        if USER_input =='c':
            while True:
                try:
                    det = input("请输入矩阵：").strip()
                    if det =='':
                        break
                    det = list(map(str, det.split(' ')))
                    det = [fractions.Fraction(item) for item in det]
                    inverse(det)
                    print('')
                    user_input = input("按回车重新运行，输入'q'退出，输入'b'回到首页：").strip().lower()
                    print(line)
                    if user_input == 'q' or user_input == 'b':
                        break
                except Exception as e:
                    print(f'Error: {e}:')
        if USER_input =='d':
            while True:
                try:
                    det1 = input("请输入左矩阵：").strip()
                    if det1 =='':
                        break
                    det1 = list(map(str, det1.split(' ')))
                    det1 = [fractions.Fraction(item) for item in det1]
                    r1 = input('行数：').strip()
                    if r1 =='':
                        break
                    c1 = input('列数：').strip()
                    if c1 =='':
                        break

                    r1 = int(r1)
                    c1 = int(c1)
                    if r1*c1!=len(det1):
                        print('矩阵或行列数有误，请重新输入\n')
                    else:
                        while True:
                            t = 0
                            times = 0
                            power = False
                            while True:
                                det2 = input("请输入右矩阵（若计算乘方则输入幂次）：").strip()
                                if det2 =='':
                                    break
                                if ' ' not in det2 and r1 == c1:
                                    if float(det2)%1 != 0:
                                        print('请输入整数\n')
                                    else:
                                        times = int(det2)
                                        det2 = det1
                                        power = True
                                        break
                                else:
                                    break
                            if det2 =='':
                                break
                            if ' ' not in det2 and r1 == c1:
                                r2 = r1
                                c2 = c1
                            else:
                                try:
                                    det2 = list(map(str, det2.split(' ')))
                                    det2 = [fractions.Fraction(item) for item in det2]
                                except Exception as e:
                                    print(f'无法计算：{e}')
                                r2 = input('行数：').strip()
                                if r2 == '':
                                    break
                                c2 = input('列数：').strip()
                                if c2 =='':
                                    break
                                try:
                                    r2 = int(r2)
                                    c2 = int(c2)
                                except Exception as e:
                                    print(f'请输入数字：{e}')
                            if r2*c2!=len(det2):
                                print('矩阵或行列数有误，请重新输入\n')
                            else:
                                if c1!=r2:
                                    print('该组矩阵不能相乘')
                                else:
                                    Product_m ={}
                                    matrix_product = None
                                    if power:
                                        l_d = r1 * c1
                                        E = [0]*l_d
                                        for i in range(r1):
                                            E[i+r1*i] = 1
                                        if times > 0 or times < 0:
                                            if times == 1:
                                                matrix_product = matrix_multiply(det1,E,t)
                                            for i in range(abs(times)-1):
                                                t += 1
                                                matrix_product = matrix_multiply(det1,det2,t,power)
                                                det2 = matrix_product(Product_m,r1)
                                            if times < 0:
                                                    inverse(det2)
                                        else:
                                            matrix_product = matrix_multiply(E,E,t-1,True)
                                    else:
                                        matrix_product = matrix_multiply(det1, det2, t)
                                    print('')
                                    while True:
                                        user_input = input("\n按回车重新运行矩阵乘法，输入'q'退出，输入'B'回到首页\n（或输入字母继续用矩积求解：a.行列式；b.伴随矩阵；c.逆矩阵；d.矩阵乘法；F.基础解系）：").strip()
                                        if user_input == 'a':
                                            try:
                                                determinant_gauss(matrix_product(Product_m,r1))
                                            except Exception as e:
                                                print(f'无法计算：{e}')
                                            user_input = input("\n计算完毕，输入'f'返回上一页，输入'q'退出，输入'B'回到首页：")
                                        if user_input == 'b':
                                            try:
                                                adjoint(matrix_product(Product_m,r1))
                                            except Exception as e:
                                                print(f'无法计算：{e}')
                                            user_input = input("\n计算完毕，输入'f'返回上一页，输入'q'退出，输入'B'回到首页：")
                                        if user_input == 'c':
                                            try:
                                                inverse(matrix_product(Product_m,r1))
                                            except Exception as e:
                                                print(f'无法计算：{e}')
                                            user_input = input("\n计算完毕，输入'f'返回上一页，输入'q'退出，输入'B'回到首页：")
                                        if user_input == 'F':
                                            try:
                                                fss(matrix_product(Product_m,r1))
                                            except Exception as e:
                                                print(f'无法计算：{e}')
                                            user_input = input("\n计算完毕，输入'f'返回上一页，输入'q'退出，输入'B'回到首页：")
                                        if user_input == 'd':
                                            det2 = matrix_product(Product_m,r1)
                                            print(line)
                                            while True:
                                                det1 = input('请输入新矩阵：').strip()
                                                if det1 =='':
                                                    break
                                                try:
                                                    det1 = list(map(str, det1.split(' ')))
                                                    det1 = [fractions.Fraction(item) for item in det1]
                                                except Exception as e:
                                                    print(f'无法计算：{e}')
                                                if det1:
                                                    position = input("请输入新矩阵位置（输入‘左/右’或'l/r'）：")
                                                    if position =='':
                                                        break
                                                    try:
                                                        if position == '左' or position == 'l':
                                                            r2 = r1
                                                            r1 = input('行数：').strip()
                                                            c1 = input('列数：').strip()
                                                            try:
                                                                r1 = int(r1)
                                                                c1 = int(c1)
                                                            except Exception as e:
                                                                print(f'请输入数字：{e}')
                                                                break
                                                            if r1 * c1 != len(det1):
                                                                print('矩阵或行列数有误，请重新输入\n')
                                                            elif c1!=r2:
                                                                user_input=''
                                                                print('\n无法计算')
                                                                break
                                                            else:
                                                                matrix_multiply(det1,det2,t)
                                                                break
                                                        if position == '右' or position == 'r':
                                                            c1 = c2
                                                            r2 = input('行数：').strip()
                                                            c2 = input('列数：').strip()
                                                            try:
                                                                r2 = int(r2)
                                                                c2 = int(c2)
                                                            except Exception as e:
                                                                print(f'请输入数字：{e}')
                                                                break
                                                            if r2 * c2 != len(det1):
                                                                print('矩阵或行列数有误，请重新输入\n')
                                                            elif c1 != r2:
                                                                user_input = ''
                                                                print('\n无法计算')
                                                                break
                                                            else:
                                                                matrix_multiply(det2,det1,t)
                                                                break
                                                    except Exception as e:
                                                        print(f'无法计算：{e}')
                                                        break
                                        if user_input == 'q' or user_input == 'B' or user_input == '':
                                            break
                                if user_input == 'q' or user_input == 'B' or user_input=='':
                                    break
                            if r2 == '' or c2 == '':
                                break
                        if det2 == '':
                            break
                        if user_input == 'q' or user_input == 'B':
                            break
                except Exception as e:
                    print(f'Error: {e}')
        if USER_input == 'e':
            while True:
                try:
                    det = input("请输入矩阵（未知数为x）：").strip()
                    if det == '':
                        break
                    det = list(map(str, det.split(' ')))
                    det_val = input('请输入行列式值：').strip()
                    if det_val == '':
                        det_val = 0
                    while True:
                        try:
                            det_val = float(det_val)
                            break
                        except Exception as e:
                            print(f'Error: {e}\n')
                            break
                    unknown_det(det,det_val)
                    print('')
                    user_input = input("按回车重新运行，输入'q'退出，输入'b'回到首页：").strip().lower()
                    print(line)
                    if user_input == 'q' or user_input == 'b':
                        break
                except Exception as e:
                    print(f'Error: {e}')
        if USER_input == 'f':
            while True:
                try:
                    det = input("请输入矩阵：").strip()
                    if det == '':
                        break
                    det = list(map(str, det.split(' ')))
                    det = [fractions.Fraction(item) for item in det]
                    col_num = input('若矩阵非方阵，请输入未知数个数（即列数） : ').strip()
                    if col_num == '':
                        fss(det)
                    else:
                        col_num = int(col_num)
                        complement = col_num - len(det) / col_num
                        det += [0] * col_num * int(complement)
                        fss(det)
                    print()
                    user_input = input("按回车重新运行，输入'q'退出，输入'b'回到首页：").strip().lower()
                    print(line)
                    if user_input == 'q' or user_input == 'b':
                        break
                except Exception as e:
                    print(f'Error: {e}')
        if USER_input == 'g':
            while True:
                try:
                    det = input("请输入矩阵：").strip()
                    if det == '':
                        break
                    det = list(map(str, det.split(' ')))
                    det_unknown = det.copy()
                    l = len(det)**(1/2)
                    if l%1 != 0:
                        print('请输入方阵!')
                    else:
                        l = int(l)
                        det_val = 0
                        for i in range(l):
                            det_unknown[i+i*l]  += '-x'
                        unknown_det(det_unknown,det_val)
                        print('\n以上为粗略特征值')
                        while True:
                            lam = input('请输入校正后的特征值：')
                            if  lam =='':
                                break
                            try:
                                det =  list(map(float, det))
                                lam = list(map(float, lam.split(' ')))
                                for i in range(len(lam)):
                                    det_lambda = det.copy()
                                    print(f'\n对特征值{i+1}（{lam[i]}）：\n')
                                    for j in range(l):
                                        det_lambda[j+j*l] -= lam[i]
                                    fss(det_lambda)
                            except Exception as e:
                                print(f'Error: {e}')
                            user_input = input("\n按回车重新输入特征值，输入'q'退出，输入'b'回到首页：").strip().lower()
                            if user_input == 'q' or user_input == 'b':
                                break
                        print(line)
                        if user_input == 'q' or user_input == 'b':
                            break
                except Exception as e:
                    print(f'Error: {e}')

        if user_input=='q':
            break
