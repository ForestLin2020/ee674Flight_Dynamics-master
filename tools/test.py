#!/usr/bin/python

#total = 0; # This is global variable.
# Function definition is here
def sum( arg1, arg2 ):
   # Add both the parameters and return them."
   total = arg1 + arg2; # Here total is local variable.
   print ("Inside the function local total : ", total)
   print(inertia()[2])
   return total;

# Now you can call sum function
total = sum( 10, 20 );
print ("Outside the function global total :",total)




def inertia():
    L  = 0000
    L1 = 111
    L2 = 222
    L3 = 333
    L4 = 444
    L5 = 555
    L6 = 555
    L7 = 777
    L8 = 888
    intr = [L, L1, L2, L3, L4, L5, L6, L7, L8]
    return intr

