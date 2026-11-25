from sympy import *
import copy

def prime_factor_decomp(gen_list):
    dict_of_factors = {}
    size = len(gen_list)
    dim = 0
    for i,gen in enumerate(gen_list):
        factors = factorint(gen)
        for factor in factors:
            if factor not in dict_of_factors:
                dict_of_factors[factor] = {}
            vec = zeros(size)[:,0]
            vec[i] = int(gen/factor)
            dim += 1
            if factors[factor] in dict_of_factors[factor]:
                dict_of_factors[factor][factors[factor]].append(vec)
            else:
                dict_of_factors[factor][factors[factor]]=  [vec]
    list_of_factors =  list(dict_of_factors.keys())
    all_factors = []
    B = zeros(size,dim)
    dim_count = 0
    index_list = []
    count = 0
    for key in sorted(list_of_factors):
        list_of_powers = dict_of_factors[key].keys()
        index_list.append(count)
        for power in sorted(list_of_powers):
            for vect in dict_of_factors[key][power]:
                B[: , dim_count] = vect
                dim_count += 1
                all_factors.append(key**power)
                count+=1
    index_list.append(count)
    return [all_factors , B,index_list]

#def diagonalize_primal_form(gen_list , Q):


def almost_diagonalize_form(gen_list , Q, index_list):
    #og_Q = copy.deepcopy(Q)
    size = len(gen_list)

    B = eye(size)
    two_dim_flag = False
    for i in range(len(index_list)-1):
        if gen_list[0]%2 == 0:
            continue
        sub_Q = Q[index_list[i]:index_list[i+1],index_list[i]:index_list[i+1]]
        sub_gen = gen_list[index_list[i]:index_list[i+1]]
        sub_size = index_list[i+1]- index_list[i]
        count = 0
        sub_B = eye(sub_size)
        two_dim_flag = False
        for factor in reversed(sub_gen):
            temp_B = eye(sub_size)#-count)
            if two_dim_flag == True:
                two_dim_flag = False
                continue
            for j in range(sub_size -count - 1, 0 , -1):
                if sub_Q[j, sub_size -count - 1].denominator == sub_gen[sub_size -count - 1]:
                    if sub_gen[0]%2 == 0 and j != sub_size -count - 1:
                        temp_vec = copy.deepcopy(temp_B[: ,sub_size -count - 2])
                        temp_B[: ,sub_size -count - 2] = copy.deepcopy(temp_B[: , j])
                        temp_B[: , j] = temp_vec
                        two_dim_flag = True
                    else:
                        temp_B[j, sub_size -count - 1] = 1
                    break
            sub_B = sub_B*temp_B

            new_sub_Q= (temp_B.T*sub_Q*temp_B)[:(sub_size-count),:(sub_size-count)]
            new_sub_gen = sub_gen[:(sub_size-count)]
            if two_dim_flag == True:
                new_B = annihilator_two_dim(new_sub_gen , new_sub_Q)
            else:
                new_B = annihilator_one_dim(new_sub_gen , new_sub_Q)
            ext_new_B = diag(new_B , eye(count))
            sub_B = sub_B*ext_new_B

            count += 1

    #    temp_B = temp_B*annihilator_one_dim(sub_gen , sub_Q)
        B[index_list[i]:index_list[i+1],index_list[i]:index_list[i+1]] = sub_B

    return B


def full_process(gen_list , Q, full_info = True):
    og_Q = copy.deepcopy(Q)
    [new_list, B1 , index_list] = prime_factor_decomp(gen_list)
    Q = B1.T*Q*B1
    Q = vector_mod_Z(Q)
    B = almost_diagonalize_form(new_list , Q, index_list)
    Q = B.T*Q*B
    return (new_list , Q , B)


def annihilator_two_dim(gen_list , sub_Q):
    #print(gen_list)
    Q = sub_Q
    size = len(gen_list)
    denom = gen_list[-1]
    rat = Rational(1,denom)
    if 2*sub_Q[-1,-1].denominator == denom and 2*sub_Q[-2,-2].denominator == denom:
        tar_Q = diag(Q[:-2,:-2] , Matrix([[2*rat,rat],[rat,2*rat]]))
    else:
        tar_Q = diag(Q[:-2,:-2] , Matrix([[0,rat],[rat,0]]))

#
#    return B
    return


def annihilator_one_dim(gen_list , sub_Q):
    #print(gen_list)
    Q = sub_Q
    size = len(gen_list)
    denom = gen_list[-1]
    num= 1
    for i in range(size):
        num = lcm(num , (Q**-1)[i , size-1].denominator)

    vec = zeros(size)[:,0]
    vec[-1] = Rational(num, denom)
    new_basis = Q**-1*vec
    B = eye(size)
    B[:,-1] = new_basis
    return B

def vector_mod_Z(mat): #gives the fractional part of the components of the vector in matrix form
    for i,comp in enumerate(mat):
        mat[i] = comp - floor(comp)
    return mat
