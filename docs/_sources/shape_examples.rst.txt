Complex structure with PoreKit() 
================================
.. code:: ipython3

    import porems as pms

Shape 1 (Single Pore)
---------------------

Create a pore with a diameter of 5 nm, a length of 10 nm, and a
reservoir length of 5 nm. The silanol on the inner surface is 6.06
:math:`\mathrm{\mu mol m^{-2}}`. The outer surface is completely
occupied by TMS molecules with 6.06 :math:`\mathrm{\mu mol m^{-2}}`.

.. code:: ipython3

    # Set pore kit
    pore = pms.PoreKit()
    
    # Set Beta-Cristobalit block witgh (x y z) (11 11 10)
    pore.structure(pms.BetaCristobalit().generate([11, 11, 10], "z"))
    pore.build()
    
    # Set a reservoir for 5 nm and an outer silanol concentration of 6.06 mumol/m2
    pore.exterior(5, hydro= 6.06)
    
    # Drill a cyclinder in the SiO2 block with a diameter of 5 and a lenght of nm with inner silanol concentration of 6.06 mumol/m2
    pore.add_shape(pore.shape_cylinder(5, 10, [5.5,5.5,5]), section={"x": [], "y": [], "z": [0,  10]}, hydro= 6.06)
    pore.prepare()
    
    # Replace silanol on the outer surface with TMS 
    pore.attach(pms.gen.tms(), mount=0, axis=[1, 2], amount=6.06, site_type="ex", inp="molar", scale=0.5)
    
    # Create the structure
    pore.finalize()
    
    # Store the structure
    pore.store("pores/shape1/")
    

.. figure::  /pics/shapes/shape1.pdf
      :align: center
      :width: 70%
      :name: fig1

Shape 2 (Parallel Pore 5 nm)
----------------------------

| Create 2 parallel pores with a diameter of 5 nm, a length of 10 nm and a reservoir length of 5 nm. The silanol on the inner surface is 6.06 :math:`\mathrm{\mu mol m^{-2}}`. The outer surface is completely occupied by TMS molecules with 6.06 :math:`\mathrm{\mu mol m^{-2}}`.

.. code:: ipython3

    # Set pore kit
    pore = pms.PoreKit()
    
    # Set Beta-Cristobalit block with (x y z) (11 11 10)
    pore.structure(pms.BetaCristobalit().generate([11, 11, 10], "z"))
    pore.build()
    
    # Set a reservoir for 5 nm and an outer silanol concentration of 6.06 mumol/m2
    pore.exterior(5, hydro= 6.06)
    
    # Drill two cyclinder in the SiO2 block with a diameter of 5 and a lenght of nm with inner silanol concentration of 6.06 mumol/m2
    pore.add_shape(pore.shape_cylinder(5, 10, [3, 3, 5]), hydro= 6.06)
    pore.add_shape(pore.shape_cylinder(5, 10, [8, 8, 5]), hydro= 6.06)
    pore.prepare()
    
    # Replace silanol on the outer surface with TMS 
    pore.attach(pms.gen.tms(), mount=0, axis=[1, 2], amount=6.06, site_type="ex", inp="molar", scale=0.5)
    
    # Create the structure
    pore.finalize()
    
    # Store the structure
    pore.store("pores/shape2/")


.. figure::  /pics/shapes/shape2.pdf
      :align: center
      :width: 70%
      :name: fig1


Shape 3 (Pore with Connection)
------------------------------

| Generate a 2 parallel pores with a diameter of 5 nm, a length of 10 nm connected with a 2.5 nm pore and a reservoir length of 5 nm. The silanol on the inner surface is 6.06 :math:`\mathrm{\mu mol m^{-2}}`. The outer surface is completely occupied by TMS molecules with 6.06 :math:`\mathrm{\mu mol m^{-2}}`.

.. code:: ipython3

    # Set pore kit
    pore = pms.PoreKit()
    
    # Set Beta-Cristobalit block witgh (x y z) (11 11 10)
    pore.structure(pms.BetaCristobalit().generate([11, 11, 10], "z"))
    pore.build()
    
    # Set a reservoir for 5 nm and an outer silanol concentration of 6.06 mumol/m2
    pore.exterior(5, hydro= 6.06)
    
    # Connect these pores
    pore.add_shape(pore.shape_cylinder(2.5, 3, [5.5, 5.5, 5], central = [1,1,0]))
    
    # Drill two cyclinder in the SiO2 block with a diameter of 5 and a lenght of nm with inner silanol concentration of 6.06 mumol/m2
    pore.add_shape(pore.shape_cylinder(5, 10, [8, 8, 5]), hydro= 6.06)
    pore.add_shape(pore.shape_cylinder(5, 10, [3, 3, 5]), hydro= 6.06)
    
    # Connect these pores
    #pore.add_shape(pore.shape_cylinder(2.5, 3, [5.5, 5.5, 5], central = [1,1,0]), hydro= 6.06)
    pore.prepare()
    
    # Replace silanol on the outer surface with TMS 
    pore.attach(pms.gen.tms(), mount=0, axis=[1, 2], amount=6.06, site_type="ex", inp="molar", scale=0.5)
    
    # Create the structure
    pore.finalize()
    
    # Store the structure
    pore.store("pores/shape3/")

.. note::
    If there are intersecting shapes, the free Si sites are assigned to the first defined shape. 
    Keep this in mind if you want to functionalize one of these shapes. 
    
.. figure::  /pics/shapes/shape3.pdf
      :align: center
      :width: 70%
      :name: fig1

Shape 4 (Pore with different inner surfaces)
--------------------------------------------

| Generate a 3 pores with a diameter of 5 nm directly adjacent to each other in z-direction. The silanol on the inner surface is 6.06 :math:`\mathrm{\mu mol m^{-2}}` in the second pore (shape_01). Shape_00 and shape_02 are occupied by TMS molecules with 6.06 :math:`\mathrm{\mu mol m^{-2}}`. The outer surface is completely occupied by TMS molecules with 6.06 :math:`\mathrm{\mu mol m^{-2}}`.

.. code:: ipython3

    # Set pore kit
    pore = pms.PoreKit()
    
    # Set Beta-Cristobalit block with (x y z) (11 11 25)
    pore.structure(pms.BetaCristobalit().generate([11, 11, 25], "z"))
    pore.build()
    
    # Set a reservoir for 5 nm and an outer silanol concentration of 6.06 mumol/m2
    pore.exterior(5, hydro= 6.06)
    
    # Drill three cyclinder with the same central in x,y direction in the SiO2 block with a diameter of 5 nm and a inner silanol concentration of 6.06 mumol/m2
    pore.add_shape(pore.shape_cylinder(5.0, 10, [5.5,5.5,5]), hydro= 6.06)   #shape_00
    pore.add_shape(pore.shape_cylinder(5.0, 5, [5.5,5.5,12.5]), hydro= 6.06) #shape_01
    pore.add_shape(pore.shape_cylinder(5.0, 10, [5.5,5.5,20]), hydro= 6.06)  #shape_02
    pore.prepare()
    
    # Replace silanol on the outer surface with TMS 
    pore.attach(pms.gen.tms(), mount=0, axis=[1, 2], amount=6.06, site_type="ex", inp="molar", scale=0.5)
    
    # Replace silanol on the inner surface with TMS for shape 0 and shape 2
    pore.attach(pms.gen.tms(), mount=0, axis=[1, 2], amount=6.06, shape="shape_00", site_type="in", inp="molar", scale=0.5)
    pore.attach(pms.gen.tms(), mount=0, axis=[1, 2], amount=6.06, shape="shape_02", site_type="in", inp="molar", scale=0.5)
    
    
    # Create the structure
    pore.finalize()
    
    # Store the structure
    pore.store("pores/shape4/")


.. figure::  /pics/shapes/shape4.pdf
      :align: center
      :width: 100%
      :name: fig1


Shape 5 (Pore with different inner surfaces)
--------------------------------------------

| Generate a 3 pores with a diameter of 5, 2.5 and 5 nm directly adjacent to each other in z-direction connected with a cone shape. The outer surface is completely occupied by TMS molecules with 6.06 :math:`\mathrm{\mu mol m^{-2}}`.

.. code:: ipython3

    # Set pore kit
    pore = pms.PoreKit()
    
    # Set Beta-Cristobalit block with (x y z) (11 11 30)
    pore.structure(pms.BetaCristobalit().generate([11, 11, 30], "z"))
    pore.build()
    
    # Set a reservoir for 5 nm and an outer silanol concentration of 6.06 mumol/m2
    pore.exterior(5, hydro= 6.06)
    
    # Drill three cyclinder with the same central in x,y direction in the SiO2 block with a diameter of 5 nm and a inner silanol concentration of 6.06 mumol/m2
    pore.add_shape(pore.shape_cylinder(5, 8, [5.5,5.5, 4]), hydro=6.06)
    pore.add_shape(pore.shape_cone(4, 3, 2,  [5.5,5.5, 9]), hydro=6.06)
    pore.add_shape(pore.shape_cylinder(2.5, 10, [5.5,5.5, 15]),  hydro=6.06)
    pore.add_shape(pore.shape_cone(3.4, 4.1, 2,  [5.5,5.5, 21]),  hydro=6.06)
    pore.add_shape(pore.shape_cylinder(5, 8, [5.5,5.5, 26]), hydro=6.06)
    pore.prepare()
    
    # Replace silanol on the outer surface with TMS 
    pore.attach(pms.gen.tms(), mount=0, axis=[1, 2], amount=6.06, site_type="ex", inp="molar", scale=0.5)
    
    # Create the structure
    pore.finalize()
    
    # Store the structure
    pore.store("pores/shape5/")


.. figure::  /pics/shapes/shape5.pdf
      :align: center
      :width: 100%
      :name: fig1

