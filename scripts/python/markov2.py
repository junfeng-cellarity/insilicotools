#!/usr/bin/env python
import random
A_count = 0
B_count = 0
C_count = 0
for i in range(1000000):
    #A shoot C
    prob = random.randint(1,100)
    if prob <= 30: #C die
        while True:
            prob = random.randint(1,100)
            if prob <= 50: #A die
                B_count +=1
                break
            else:
                prob = random.randint(1,100)
                if prob > 50: #B die
                    A_count += 1
                    break
    else:#B die
        prob = random.randint(1,100)
        if prob > 50: #C die
            while True:
                prob = random.randint(1,100)
                if prob <= 30: #B die
                    A_count +=1
                    break
                else:
                    prob = random.randint(1,100)
                    if prob > 50: #A die
                        B_count += 1
                        break
        else:#B die
            prob = random.randint(1,100)#A shoot C
            if prob <= 30:#C die
                A_count += 1
            else:
                C_count += 1
print A_count,B_count,C_count


        # C die

