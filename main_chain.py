# In this varible we keep the final name of the molecule
name_of_the_molecule = ""

# This dictionary describe the main chain of the molecule
chemical_prefixes = {1: "meth", 2: "eth", 3: "prop", 4: "but", 5: "pent",
                    6: "hex", 7: "hept", 8: "oct", 9: "non", 10: "dec",
                    11: "undec", 12: "docec", 13: "tridec", 14: "tetradec", 15: "pentadec",
                    16: "hexadec", 17: "heptadec", 18: "octadec", 19: "nonadec", 20: "icos",
                    21: "henicos", 22: "docos", 23: "tricos", 24: "tetracos", 25: "pentacos",
                    26: "hexacos", 27: "heptacos", 28: "octacos", 29: "nonacos", 30: "triacont",
                    31: "hentriacont", 32: "dotriacont", 33: "tritriacont", 34: "tetratriacont", 35: "pentatriacont",
                    36: "hexatriacont", 37: "heptatriacont", 38: "octatriacont", 39: "nonatriacont", 40: "tetracont",
                    41: "hentetracont", 42: "dotetracont", 43: "tritetracont", 44: "tetratetracont", 45: "pentatetracont",
                    46: "hexatetracont", 47: "heptatetracont", 48: "octatetracont", 49: "nonatetracont", 50: "pentacont",
                    51: "henpentacont", 52: "dopentacont", 53: "tripentacont", 54: "tetrapentacont", 55: "pentapentacont",
                    56: "hexapentacont", 57: "heptapentacont", 58: "octapentacont", 59: "nonapentacont",
                    60: "hexacont", 61: "henhexacont", 62: "dohexacont", 63: "trihexacont", 64: "tetrahexacont", 65: "pentahexacont",
                    66: "hexahexacont", 67: "heptahexacont", 68: "octahexacont", 69: "nonahexacont", 70: "heptacont",
                    71: "henheptacont", 72: "doheptacont", 73: "triheptacont", 74: "tetraheptacont", 75: "pentaheptacont",
                    76: "hexaheptacont", 77: "heptaheptacont", 78: "octaheptacont", 79: "nonaheptacont", 80: "octacont",
                    81: "henoctacont", 82: "dooctacont", 83: "trioctacont", 84: "tetraoctacont", 85: "pentaoctacont",
                    86: "hexaoctacont", 87: "heptaoctacont", 88: "octaoctacont", 89: "nonaoctacont", 90: "enneacont",
                    91: "henenneacont", 92: "doenneacont", 93: "trienneacont", 94: "tetraenneacont", 95: "pentaenneacont",
                    96: "hexaenneacont", 97: "heptaenneacont", 98: "octaenneacont", 99: "nonaenneacont", 100: "hect"}

def iupac(name):

    # The function calls the string in the dictionary chemical_prefixes
    def prefix(carbons):
        global chemical_prefixes
        if carbons in chemical_prefixes:
            return chemical_prefixes[carbons]
    
    # This fucntion inverse the string to read the smiles 
    def inverse(string):
        string = string.replace("(","~")
        string = string.replace(")","(")
        string = string.replace("~",")")
        return string[::-1]
    
    # The function reduce the string to numbers for read the SMILES and look patterns CC(C)CCC -> 2(1)3
    def encode_c_string(s):
        result = []
        current_count = 0
        s = s + "~"
        for char in s:
            if char == 'C':
                current_count += 1
            elif char == "(" :
                if current_count > 0:
                    result.append(current_count)
                    result.append("(")
                    current_count = 0
                else:
                    result.append("(")
                    current_count = 0
            elif char == ")":
                if current_count > 0:
                    result.append(current_count) 
                    result.append(")")
                    current_count = 0
                else:
                    result.append(current_count)
                    result.append(")")
                    current_count = 0
            elif char == "~":
                if current_count > 0:
                    result.append(current_count)
                    current_count = 0
        return result
    
    # This function name the alkyl 
    def mainchain(chain):
        carbons = 0
        chain.append("~")
        main_chain = 0

        # This case is just when the chain don't have branches 
        if len(chain) == 2:
            carbons += chain[0]
        
        # The chain is not lineal
        else:

            # With this while the for loop follow until the list chain is one
            while len(chain) > 1:

                # The for loop read the chain from left to right 
                for char in range(len(chain)):

                    # Case a(b)c~ where a, b, c are int numbers
                    if chain[char] == ")" and len(chain) == 6:
                        # c > b or c > a where c is not 1
                        if chain[char+1] > chain[char-3] or chain[char+1] > chain[char-1]:
                            # a > b
                            if chain[char-3] > chain[char-1]:
                                carbons =+ chain[char-3] + chain[char+1]
                                chain[char+1] = carbons
                                for i in range(4):
                                    chain.pop(char-3)
                                break
                            # b > a or b = a
                            else:
                                carbons =+ chain[char-1]  + chain[char+1]
                                chain[char+1] = carbons
                                for i in range(4):
                                    chain.pop(char-3)
                                break
                        # a(b)1~ b > c or a > c, c is 1 
                        elif chain[char+1] == 1:
                            carbons =+ chain[char-3] + chain[char-1] + 1 
                            chain[char+1] = carbons
                            for i in range(4):
                                chain.pop(char-3)
                            break
                        # b > c or a > c, where c is not 1
                        else:
                            if chain[char-3] > chain[char-1]:
                                carbons =+ chain[char-3] + chain[char-1] + 1                      
                                chain[char+1] = carbons
                                for i in range(4):
                                    chain.pop(char-3)
                                break
                            else:
                                carbons =+ chain[char-3] + chain[char-1] + 1                      
                                chain[char+1] = carbons
                                for i in range(4):
                                    chain.pop(char-3)
                                break
                    
                    # Case (a)(b) where a and b are int numbers
                    elif chain[char] == ")" and chain[char+1] == "(" and isinstance(chain[char+2], int) and chain[char+3] == ")":
                        if (chain[char-1] + chain[char+2]) > main_chain:
                            main_chain = chain[char-1] + chain[char+2] + 1
                            break
                        # c((a)(b)d)e                          
                        elif chain[char-3] == "(":
                            # a and b are any int numbers
                            if chain[char-1] > chain[char+2]:
                                for i in range(3):
                                    chain.pop(char+1)
                                break
                            else:
                                chain[char-1] = chain[char+2]
                                for i in range(3):
                                    chain.pop(char+1)
                                break
                        # c(a)(b)d                    
                        elif isinstance(chain[char-3], int):
                            # 1(a)(b)d, where c is 1
                            if chain[char-3] == 1:
                                # a > b
                                if chain[char-1] > chain[char+2]:
                                    a = chain[char-3]
                                    chain[char-3] = chain[char-1]
                                    chain[char-1] = a
                                    count = 0
                                    break
                                # 1(1)(b)d, where c,a  is 1
                                elif chain[char-3] == 1 and chain[char-1] == 1:
                                    for i in range(3):
                                        chain.pop(char+1)
                                    break
                                # 1(a)(1)d, where c,b  is 1
                                elif chain[char-3] == 1 and chain[char+2] == 1:
                                    for i in range(3):
                                        chain.pop(char-3)
                                    break
                                # b > a or b = a
                                else:
                                    a = chain[char-3]
                                    chain[char-3] = chain[char+2] 
                                    chain[char+2] = a
                                    break
                            
                            # c(a)(b)d, where c is any number                                                             
                            else:
                                # a > b
                                if chain[char-1] > chain[char+2]:
                                    # a > c
                                    if chain[char-1] > chain[char-3]:                                    
                                        for i in range(3):
                                            chain.pop(char+1)
                                        break
                                    # c > a or c = a                                    
                                    else:
                                        for i in range(3):
                                            chain.pop(char+1)
                                        break
                                # b > a or b = a
                                else:
                                    # b > c
                                    if chain[char+2] > chain[char-3]:
                                        chain[char-1] = chain[char+2]
                                        for i in range(3):
                                            chain.pop(char+1)
                                        break
                                    # c > b or c = b                                    
                                    else:
                                        chain[char-1] = chain[char+2]                                        
                                        for i in range(3):
                                            chain.pop(char+1)
                                        break
                    
                    # Case a(b(c)d)e where a, b, c, d, e are int numbers
                    elif chain[char] == ")" and isinstance(chain[char+1], int) and len(chain) != 6 and isinstance(chain[char-3], int):
                        #print(chain)
                        if (chain[char-1] + chain[char-3]) > main_chain:
                            main_chain = chain[char-1] + chain[char-3] + 1
                            break
                        # c > a 
                        elif chain[char-1] > chain[char-3]:
                            carbons =+ chain[char-1] + chain[char+1]
                            chain[char+1] = carbons
                            for i in range(4):
                                chain.pop(char-3)
                            break
                        # b > c or b = c
                        else:
                            carbons =+ chain[char-3] + chain[char+1]
                            chain[char+1] = carbons
                            for i in range(4):
                                chain.pop(char-3)
                            break
                    
                    # Case ((b)c) or ~(b)c
                    elif chain[char] == ")" and isinstance(chain[char+1], int) and not isinstance(chain[char-3], int):
                        chain[char+1] = chain[char-1] + chain[char+1]
                        for i in range(3):
                            chain.pop(char-2)
                        break
                    
                    # Case a~
                    elif chain[char] == "~" and len(chain) == 2:
                        if main_chain > chain[char-1]:
                            chain[char-1] = main_chain
                            carbons = main_chain
                        else:
                            chain.pop(-1)
                    
                    # Case the SMILES is not start with a C example: (CC)CC
                    elif len(chain) == 5:
                        return print("Incorrect SMILES")
                    
                    # Case No previous case fits
                    else:
                        pass
        return carbons  

    #Create a copy of the original SMILES
    chain = name

    global name_of_the_molecule
    if ("=" not in chain and "#" not in chain
        and "\\" not in chain and "/" not in chain
        and "@" not in chain and "O" not in chain
        and "N" not in chain and "S" not in chain
        and "F" not in chain and "Cl" not in chain):
        #Give the name for alcanes with this if
        chain.upper()
        name_of_the_molecule += prefix(mainchain(encode_c_string(inverse(chain))))
        name_of_the_molecule += "ane"
    
    print(name_of_the_molecule)