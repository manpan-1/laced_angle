from part import *
from material import *
from section import *
import abq_toolset as xtr
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import sys
import os
from shutil import copyfile
import odbAccess
import EN_tools as en

##############INPUT DATA ####################

# Cross section b / (epsilon * t) according to EC3-1-1 table 5.1 sheet 3 for angles
p_class = 15.

# Total width of the leg (from the intersection of the midlines to the edge)
l_leg = 100.

# flexural slenderness
lambda_flex = 0.4

# Lacing density. Lacing dist over leg width ratio
lace_over_leg = 1.1

# Imperfection amplitude for local, flexural and torsional modes
loc_imp = 200.
flex_imp = 200.
tor_imp = 200.

# Yield stress
fy_steel = 690.

try:
    l_leg = float(sys.argv[-4])
except:
    pass

try:
    p_class = float(sys.argv[-3])
except:
    pass

try:
    lambda_flex = float(sys.argv[-2]) / 100
except:
    pass

try:
    lace_over_leg = float(sys.argv[-1]) / 10
except:
    pass


# MATERIAL
# Modulus of elasticisy, poissons ratio, yield stress and epsilon
# (Yield stress used for design calculations, not for material properties of the model)
v_poi = 0.3
E_steel = 210000.
G_steel = E_steel / (2 * (1 + v_poi))
epsilon = sqrt(235 / fy_steel)

# BASIC GEOMETRY
# Shell thickness
t_shell = l_leg / (p_class * epsilon)

# The straight part of the leg (l_strght) and the corner bending radius
# The bending radius is defined as a ratio to the shell thickness
r_over_t = 3.
r_bend = r_over_t * t_shell
l_strght = l_leg - r_bend


# Create a string which will be used as an identifier for the current model. Used for directory and filename.
IDstring = str(int(l_leg))+'-'+\
           str(int(p_class))+'-'+\
           str(int(lambda_flex * 100))+'-'+\
           str(int(lace_over_leg * 10))

# Make a new subdirectory for the current session
os.mkdir(IDstring)

# Copy necessary files to the new directory
copyfile('abq_toolset.py', './'+IDstring+'/abq_toolset.py')
copyfile('laced_angle.py', './'+IDstring+'/laced_angle.py')
copyfile('EN_tools.py', './'+IDstring+'/EN_tools.py')
#copyfile('GN_Riks_killer.f', './'+IDstring+'/GN_Riks_killer.f')

# Change the working directory
os.chdir('./'+IDstring)


# CROSS SECTION PROPERTIES

# CORNER
# Outer radius
r_max = r_bend + (t_shell / 2)

# Inner radius
r_min = r_bend - (t_shell / 2)

# AREA
# Area of leg in Y-axis
A1 = l_strght * t_shell

# Area of leg in X-axis
A2 = A1

# Area of large quarter circle
A3 = r_max ** 2 * pi / 4

# Area of small quarter circle
A4 = r_min ** 2 * pi / 4

# Total area
A_tot = A1 + A2 + A3 - A4

# CENTRE OF GRAVITY
# Leg in Y-axis
yg1 = r_bend + l_strght / 2

# Leg in x-axis
yg2 = t_shell / 2

# Large quarter circle
yg3 = r_bend * (1 - 4 / (3 * pi))

# Small quarter circle
yg4 = t_shell + (r_bend - t_shell) * (1 - 4 / (3 * pi))

# Center of gravity from Origo
ytp = (A1 * yg1 + A2 * yg2 + A3 * yg3 - A4 * yg4) / A_tot

# MOMENT OF INERTIA
# Moment of intertia part 1 (leg in Y-axis)
Iz1 = t_shell * l_strght ** 3 / 12 + l_strght * t_shell * (yg1 - ytp) ** 2

# Moment of intertia part 2 (leg in X-axis)
Iz2 = l_strght * t_shell ** 3 / 12 + l_strght * t_shell * (yg2 - ytp) ** 2

# Moment of intertia part 3 (large quarter circle)
Iz3 = r_bend ** 4 * (pi / 16 - 4 / (9 * pi)) + pi * r_bend ** 2 / 4 * (yg3 - ytp) ** 2

# Moment of intertia part 3 (small quarter circle)
Iz4 = (r_bend - t_shell) ** 4 * (pi / 16 - 4 /(9 * pi)) + pi * (r_bend - t_shell) ** 2 / 4*(yg4 - ytp) ** 2

# Moments of intertia
I_z = Iz1 + Iz2 + Iz3 - Iz4
I_y = I_z

# Product of intertia
Iyz = A1 * (yg1 - ytp) * (yg2 - ytp) * 2 + A3 * (yg3 - ytp) ** 2 - A4 * (yg4 - ytp) ** 2

# Principal moments and X-Y direction about centroid
I1 = I_z + abs(Iyz)
I2 = I_z - abs(Iyz)

# Bending stiffness E*I
E_I1 = E_steel * I1
E_I2 = E_steel * I2

# COLUMN LENGTH
l_tot = lambda_flex * pi * sqrt(E_I2 / (A_tot * fy_steel))

# LACING
# Approximate longitudinal spacing of lacing attachment points
l_lace_approx = lace_over_leg * l_leg

# Calculate the number of laces that best fit to the requested ratio
n_laces = int(round(l_tot / l_lace_approx))

# adjust l_tot
l_lace = l_lace_approx
l_tot = n_laces * l_lace

# Final lacing spacing
l_lace = l_tot / n_laces

# Calculate the actual lacing space over leg width ratio
real_lace_over_leg = l_lace / l_leg

# The following calculations conclude to the cross-sectional area of the lacing rods.
# It is calculated for a tubular section of slenderness = 1 and d/t = 70

# distance of two opposite points on the legs' edges
l_hypot = sqrt(2) * l_leg

# Rod length
l_wire = sqrt(l_lace ** 2 + l_hypot ** 2)

# Rod slenderness
lambda_wire = 1.

# Rod class
wire_classification = 70.

# the following two constants (constant a and b) are used to assist the calculations
constant_a = (lambda_wire * pi / l_wire) ** 2 * E_steel / fy_steel
constant_b = 2 / (70 * epsilon ** 2)

# Outer radius of the tube
r_wire = sqrt(4 * (1 - (1 - constant_b) ** 2) / (constant_a * (1 - (1 - constant_b) ** 4)))

# Area.
# For the case of no lacing, uncomment the second line to decrease the laceing cs are
A_wire = pi * r_wire ** 2 * (1 - (1 - constant_b) ** 2)
#A_wire = 1.e-8

# DESIGN RESISTANCE
# The following calculations take are for the design resistance of the compression element
# according to EC1-1-1 and EC3-1-5 and include overall and local buckling resistance.
# The calculations regard to a plain L-profile neglecting the lacing.


# Classification and Aeff
# Aeff is calculated assuming uniform compression on the sectors.
psi = 1.
kapa_sigma = 0.57 - 0.21 * psi + 0.07 * psi ** 2
p_class = l_leg / (t_shell * epsilon)
lambda_p = p_class / (28.4 * sqrt(kapa_sigma))
if lambda_p > 0.748 and int(p_class) > 15:
    rho = (lambda_p - 0.188) / lambda_p ** 2
else:
    rho = 1.

# Effective cross-section
A_eff = 2 * A1 * rho + A3 - A4

# Overall buckling.Calculation of Ncr.
# All three modes (bending over the principal axes and torsional) are taken into account.
# Formulas from Trahair.
I_torsion = 2 * (t_shell ** 3 * l_leg / 3)
I_warp = 0.

# Calculate the critical flexural-torsional load using the equivalent
# function from EN_tools module.
N_cr = en.N_cr_flex_tor(l_tot, A_eff, I_y, I_z, Iyz, I_torsion, I_warp, y_sc=ytp, z_sc=ytp)

# Design plastic resistance
N_pl_rd = fy_steel * A_eff

# Member slenderness
lambda_flex_tor = sqrt(N_pl_rd / N_cr)

# Reduction factor chi
a_imp_fact = 0.34
phi_capital = (1 + a_imp_fact * (lambda_flex_tor - 0.2) + lambda_flex_tor ** 2) / 2
chi_glob = 1 / (phi_capital + sqrt(phi_capital ** 2 - lambda_flex_tor ** 2))

# Buckling resistance
N_b_rd = chi_glob * N_pl_rd

# WRITE OUT FILE
# model information are writen in a text file
out_file = open('./model_info-'+IDstring+'.dat', 'w')
out_file.write('\n-GEOMETRIC CHARACTERISTICS\n')
out_file.write('Total leg width:..................................................... '+str(l_leg)+' [mm]\n')
out_file.write('Bending radius (midline):............................................ '+str(r_bend)+' [mm]\n')
out_file.write('Length of the leg\'s flat part:....................................... '+str(l_strght)+' [mm]\n')
out_file.write('Total column length:................................................. '+str(l_tot/1000)+' [m]\n')
out_file.write('Profile thickness:................................................... '+str(t_shell)+' [mm]\n')
out_file.write('Number of lacing bars:............................................... '+str(int(n_laces))+'\n')
out_file.write('Lacing length over leg width:........................................ '+str(lace_over_leg)+'\n')
out_file.write('\n-STRUCTURAL CHARACTERISTICS'+'\n')
out_file.write('Yield strength:...................................................... '+str(fy_steel)+' [MPa]\n')
out_file.write('Gross cross-sectional area:.......................................... '+str(A_tot)+' [mm^2]\n')
out_file.write('Cross-sectional area of lacing rods, A_rod:.......................... '+str(A_wire)+' [mm^2]\n')
out_file.write('Moment of inertia aroung the axes parallel to the legs, Iy, Iz:...... '+str(I_y)+' [mm^4]\n')
out_file.write('Max principal moment of inertia, I1:................................. '+str(I1)+' [mm^4]\n')
out_file.write('Min principal moment of inertia, I2:................................. '+str(I2)+' [mm^4]\n')
out_file.write('Cross-section classification, c/(epsilon*t):......................... '+str(p_class)+'\n')
out_file.write('Plate slenderness, lambda_p:......................................... '+str(lambda_p)+'\n')
out_file.write('Effective area reduction factor, rho:................................ '+str(rho)+'\n')
out_file.write('Effective cross-sectional area:...................................... '+str(A_eff)+' [mm^2]\n')
out_file.write('Max flexural critical load, N_cr_max:................................ '+str(N_cr_max/1000)+' [kN]\n')
out_file.write('Min flexural critical load, N_cr_min:................................ '+str(N_cr_min/1000)+' [kN]\n')
out_file.write('Torsional critical load, N_cr_tor:................................... '+str(N_cr_tor/1000)+' [kN]\n')
out_file.write('Combined torsional-flexural critical load, N_cr:..................... '+str(N_cr/1000)+' [kN]\n')
out_file.write('Flexural slenderness, lambda_flex:................................... '+str(lambda_flex)+'\n')
out_file.write('Flexural-torsional slenderness, lambda_tor_flex:..................... '+str(lambda_flex_tor)+'\n')
out_file.write('Flexural-torsional buckling reduction factor, chi:................... '+str(chi_glob)+'\n')
out_file.write('Plastic resistance, N_pl_rd:......................................... '+str(N_pl_rd/1000)+' [kN]\n')
out_file.write('Buckling resistance, N_b_rd:......................................... '+str(N_b_rd/1000)+' [kN]\n')
out_file.write('\n-MODEL IMPERFECTIONS'+'\n')
out_file.write('Flexural buckling bow imperfections:................................. l/'+str(flex_imp)+'\n')
out_file.write('Tortional buckling imperfections:.................................... l/'+str(tor_imp)+'\n')
out_file.write('Plate imperfections:................................................. b/'+str(loc_imp)+'\n')

### MODEL ###
#Create a model
mdb.models.changeKey(
    fromName = 'Model-1',
    toName = 'BracedLBeam'
    )

beamModel = mdb.models['BracedLBeam']

#Sketch: Create the cross-section: L-section
p11 = (0, l_leg)
p12 = (0, 0)
p21 = (l_leg, 0)
p22 = (0, 0)
np1 = (0, l_leg)
np2 = (l_leg, 0)

cs_sketch = beamModel.ConstrainedSketch(
    name = '__profile__',
    sheetSize = l_leg + r_bend
    )

cs_sketch.Line(
    p11,
    p12
    )

cs_sketch.VerticalConstraint(
    addUndoState = False,
    entity = cs_sketch.geometry[2]
    )
    
cs_sketch.Line(
    p21,
    p22
    )

#Sketch: Create angle to the L-section
cs_sketch.HorizontalConstraint(
    addUndoState = False,
    entity = cs_sketch.geometry[3]
    )
    
cs_sketch.PerpendicularConstraint(
    addUndoState = False,
    entity1 = cs_sketch.geometry[2],
    entity2 = cs_sketch.geometry[3]
    )
    
cs_sketch.FilletByRadius(
    curve1 = cs_sketch.geometry[2],
    curve2 = cs_sketch.geometry[3],
    nearPoint1 = np1,
    nearPoint2 = np2,
    radius =  r_bend)
    
#Creating the beam: Extrudes the cross-section
column_prt = beamModel.Part(
    dimensionality = THREE_D,
    name = 'Part',
    type = DEFORMABLE_BODY
    )
    
column_prt.BaseShellExtrude(
    depth = l_tot,
    sketch = cs_sketch
    )

#del cs_sketch

#Creating datumplanes
#Creatng datumplane every l_lace in the range 0 < ext

for datumplanes in range(1, n_laces) : 
    column_prt.DatumPlaneByPrincipalPlane(
        offset = datumplanes * l_lace,
        principalPlane = XYPLANE
        )

## Creating a partition on every datumplane    
for datum_partition in column_prt.datums.items():
    column_prt.PartitionFaceByDatumPlane(
        datumPlane = datum_partition[1],
        faces = column_prt.faces[:]
        )

# Creating the zigzag bracing for the column_prt
for wir1 in range(0, n_laces, 2) :
    column_prt.WirePolyLine(
        points = ((0.,l_leg,wir1 * l_lace),(l_leg,0.,wir1 * l_lace+l_lace),),
        meshable = ON
		)

for wir2 in range(1, n_laces, 2) :
    column_prt.WirePolyLine(points = ((l_leg,0. , wir2 * l_lace),(0. ,l_leg ,wir2 * l_lace+l_lace),), meshable = ON)    


#Defining material: Steel
beamModel.Material(name = 'steel')
beamModel.materials['steel'].Elastic(table = ((E_steel, v_poi), ))
beamModel.materials['steel'].Plastic(
    table = (
        (fy_steel, 0.0), 
        (fy_steel * 1.2, 0.2)
        )
    )

#Creating section: Shell
beamModel.HomogeneousShellSection(
    idealization = NO_IDEALIZATION, 
    integrationRule = SIMPSON,
    material = 'steel',
    name = 'shell',
    numIntPts = 5, 
    poissonDefinition = DEFAULT,
    preIntegrate = OFF,
    temperature = GRADIENT, 
    thickness = t_shell,
    thicknessField = '',
    thicknessModulus = None,
    thicknessType = UNIFORM,
    useDensity = OFF
    )
    
#Creating section: Beam
beamModel.TrussSection(area = A_wire,
    material = 'steel',
    name = 'bracing'
	)

# Create set with all faces 
all_faces = column_prt.Set(
    faces = column_prt.faces[:],
    name = 'All_Faces'
    )

#Assign section: To the shell
column_prt.SectionAssignment(
    offset = 0.0, 
    offsetField = '',
    offsetType = MIDDLE_SURFACE,
    region = all_faces,
    sectionName = 'shell', 
    thicknessAssignment = FROM_SECTION
    )

# Create a set for the wires
edge_list = []
for z_lace in range(1, 2*n_laces, 2):
	point_coord = ((l_leg*0.5, l_leg*0.5, l_lace * z_lace / 2), )	
	edge_list.append(column_prt.edges.findAt(point_coord))

# The findAt function does not seem to work. The required arguments are 2. The find_points has 11 argumenst.
all_wires = column_prt.Set(
    edges = edge_list,
	name = 'all_wires'
	)	

#Assign section: To the bracing
column_prt.SectionAssignment(
    offset = 0.0, 
    offsetField = '', 
    offsetType = MIDDLE_SURFACE, 
    region = all_wires, 
    sectionName = 'bracing', 
    thicknessAssignment = FROM_SECTION
    )

#Creating material orientation: The the shell    
column_prt.MaterialOrientation(
    additionalRotationType = ROTATION_NONE, 
    axis = AXIS_2, 
    fieldName = '', 
    localCsys = None, 
    orientationType = GLOBAL, 
    region = all_faces
    )

# Assign truss element type for the lacing wires
column_prt.setElementType(
    elemTypes = (ElemType(
        elemCode = T3D2, 
        elemLibrary = STANDARD
        ), ),
    regions = all_wires
    )

#Creating the assembly
beamModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
beamAssembly = beamModel.rootAssembly

AssemblyInstance = beamAssembly.Instance(
    dependent = OFF,
    name = 'Column Instance', 
    part = column_prt
    )

#Step: Creating the riks step
beamModel.StaticRiksStep(
    initialArcInc = 0.1,
    name = 'RIKS',
    nlgeom = ON,
    previous = 'Initial',
    maxNumInc = 25
    )

#Creating RP: reference points
cog1 = column_prt.getMassProperties()
#cog1['areaCentroid']
lst = list(cog1['areaCentroid'])
lst[2] = 0
cog1 = tuple(lst)

cog2 = column_prt.getMassProperties()
#cog2['areaCentroid']
lst = list(cog2['areaCentroid'])
lst[2] = l_tot
cog2 = tuple(lst)

beamAssembly.ReferencePoint(point = cog1)
beamAssembly.ReferencePoint(point = cog2)

# Creating a Set for the selection at RP - 1, the one at the bottom of the column
beamAssembly.Set(
    name = 'RP-1', 
    referencePoints = (
    beamAssembly.referencePoints[4], ))

# Creating a Set for selecting the edges at the base of the column
beamAssembly.Set(
    edges = AssemblyInstance.edges.getByBoundingBox( 0, 0, 0, l_leg+1, l_leg+1, 0),
    name = 'Base_edge'
    )
	
#Coupling: The RPs to the ends    
beamModel.Coupling(
    controlPoint = beamAssembly.sets['RP-1'],
    couplingType = KINEMATIC, 
    influenceRadius = WHOLE_SURFACE,
    localCsys = None,
    name = 'Constraint-1', 
    surface = beamAssembly.sets['Base_edge'],
    u1 = ON, 
    u2 = ON, 
    u3 = ON, 
    ur1 = ON, 
    ur2 = ON, 
    ur3 = ON
    )

# Creating a Set at RP - 2, , the one at the top of the column
beamAssembly.Set(
    name = 'RP-2',
    referencePoints = (beamAssembly.referencePoints[5], )
    )

# Creating a Set for selecting the edges at the head of the column
beamAssembly.Set(
    edges = AssemblyInstance.edges.getByBoundingBox( 0, 0, l_tot, l_leg+1, l_leg+1, l_tot),
    name = 'Top_edge'
    )

#Coupling: The RPs to the ends
beamModel.Coupling(
    controlPoint = beamAssembly.sets['RP-2'],
    couplingType = KINEMATIC, 
    influenceRadius = WHOLE_SURFACE,
    localCsys = None,
    name = 'Constraint-2', 
    surface = beamAssembly.sets['Top_edge'],
    u1 = ON, 
    u2 = ON, 
    u3 = ON,
    ur1 = ON, 
    ur2 = ON, 
    ur3 = ON
    )
  
#Boundary condition: Bottom : boundary condition set on the RPs
beamModel.DisplacementBC(
    amplitude = UNSET,
    buckleCase = PERTURBATION_AND_BUCKLING, 
    createStepName ='RIKS', 
    distributionType = UNIFORM, 
    fieldName = '', 
    fixed = OFF, 
    localCsys = None, 
    name='BC-1', 
    region = beamAssembly.sets['RP-1'], 
    u1 = 0.0, 
    u2 = 0.0, 
    u3 = 0.0, 
    ur1 = UNSET, 
    ur2 = UNSET, 
    ur3 = 0.0 
    )
    
#Boundary condition: Top : boundary condition set on the RPs
beamModel.DisplacementBC(
    amplitude = UNSET, 
    buckleCase = PERTURBATION_AND_BUCKLING, 
    createStepName = 'RIKS', 
    distributionType = UNIFORM, 
    fieldName = '', 
    fixed = OFF, 
    localCsys = None, 
    name = 'BC-2', 
    region = beamAssembly.sets['RP-2'], 
    u1 = 0.0, 
    u2 = 0.0, 
    u3 = UNSET, 
    ur1 = UNSET, 
    ur2 = UNSET, 
    ur3 = UNSET
    )

#Mesh: Changed type to structured mesh
beamAssembly.setMeshControls(
    regions = AssemblyInstance.faces[:], 
    technique = STRUCTURED
    )

# seed the lacing wires so that they are a single element per brace
beamAssembly.seedEdgeByNumber(
    constraint = FINER, 
    edges = AssemblyInstance.sets['all_wires'].edges[:],
    number = 1
    )

# Seed the shell
beamAssembly.seedPartInstance(
    deviationFactor = 0.1, 
    minSizeFactor = 0.1, 
    regions = (AssemblyInstance, ), 
    size = l_leg / 10
    )

# Mesh: Generating the mesh
beamAssembly.generateMesh(
    regions = (AssemblyInstance, )
    )

## LOCAL AND GLOBAL BUCKLING IMPERFECTIONS	
loc_imperfection_amp = l_lace / loc_imp
flex_imperfection_amp = l_tot / flex_imp
tor_imperfection_amp = l_tot / tor_imp

beamAssembly.makeIndependent(
    instances = (
        AssemblyInstance,
        )
    )

for j in range(len(AssemblyInstance.nodes)):
    xi = AssemblyInstance.nodes[j].coordinates[0]
    yi = AssemblyInstance.nodes[j].coordinates[1]
    zi = AssemblyInstance.nodes[j].coordinates[2]
    glob_bow = -flex_imperfection_amp * sin(pi * zi / l_tot) * sin(pi / 4)
    theta_tor = tor_imperfection_amp / (l_leg + r_bend) * sin(pi * zi / l_tot)
    dx_torsion = (xi * cos(theta_tor) - yi * sin(theta_tor) - xi)
    dy_torsion = (xi * sin(theta_tor) + yi * cos(theta_tor) - yi)
    if xi==0:
        beamAssembly.editNode(
            nodes = AssemblyInstance.nodes[j],
            offset1 = glob_bow + dx_torsion + loc_imperfection_amp * sin(pi * zi / (2 * l_lace)) * ((yi - r_bend) / (l_leg - r_bend)),
            offset2 = glob_bow + dy_torsion
            )
    elif yi==0:
        beamAssembly.editNode(
            nodes = AssemblyInstance.nodes[j],
            offset1 = glob_bow + dx_torsion,
            offset2 = glob_bow + dy_torsion + loc_imperfection_amp * cos(pi / 1. + pi * (zi / (2*l_lace))) * ((xi - r_bend) / (l_leg - r_bend))
            )
    else:
        beamAssembly.editNode(
            nodes = AssemblyInstance.nodes[j],
            offset1 = glob_bow + dx_torsion,
            offset2 = glob_bow + dy_torsion
            )

# Add concentrated force
beamModel.ConcentratedForce(
    cf3 = -N_b_rd,
    createStepName = 'RIKS',
    distributionType = UNIFORM,
    field = '',
    localCsys = None,
    name = 'compression',
    region = beamAssembly.sets['RP-2']
    )


# Field and History output requests

beamModel.historyOutputRequests.changeKey(
    fromName = 'H-Output-1',
    toName = 'load'
    )
    
beamModel.historyOutputRequests['load'].setValues(
    rebar = EXCLUDE,
    region = beamAssembly.sets['RP-1'], 
    sectionPoints = DEFAULT, variables = ('RF3', )
    )
    
beamModel.HistoryOutputRequest(
    createStepName = 'RIKS',
    name = 'disp',
    rebar = EXCLUDE,
    region = beamAssembly.sets['RP-2'],
    sectionPoints = DEFAULT,
    variables = ('U3', )
    )

beamModel.fieldOutputRequests.changeKey(
    fromName = 'F-Output-1', 
    toName = 'fields'
    )
beamModel.fieldOutputRequests['fields'].setValues(
    variables = ('S', 'MISES', 'E', 'PEEQ', 'U')
    )


# Creating the riks analysis
# Job: Creating the job
riks_job = mdb.Job(
    atTime = None,
    contactPrint = OFF,
    description = '',
    echoPrint = OFF, 
    explicitPrecision = SINGLE,
    getMemoryFromAnalysis = True,
    historyPrint = OFF, 
    memory = 90,
    memoryUnits = PERCENTAGE,
    model = 'BracedLBeam',
    modelPrint = OFF, 
    multiprocessingMode = DEFAULT,
    name = 'riks-job-'+IDstring,
    nodalOutputPrecision = SINGLE, 
    numCpus = 4,
    numDomains = 4,
    numGPUs = 0,
    queue = None,
    scratch = '',
    type = ANALYSIS, 
    userSubroutine = '',
    waitHours = 0,
    waitMinutes = 0
    )

#Job: Submiting the job    
riks_job.submit(consistencyChecking = OFF)

# Wait for completion
riks_job.waitForCompletion()


# Collect the max results and write rhem in the output file
odb_name = 'riks-job-'+IDstring
myOdb = odbAccess.openOdb(path=odb_name+'.odb')
RIKSstep = myOdb.steps['RIKS']
rp1key = RIKSstep.historyRegions.keys()[1]
ho1key = RIKSstep.historyRegions[rp1key].historyOutputs.keys()[0]
rp2key = RIKSstep.historyRegions.keys()[2]
ho2key = RIKSstep.historyRegions[rp2key].historyOutputs.keys()[0]
asskey = RIKSstep.historyRegions.keys()[0]
hoasse = RIKSstep.historyRegions[asskey].historyOutputs.keys()[-1]
load_hist = RIKSstep.historyRegions[rp1key].historyOutputs[ho1key].data
disp_hist = RIKSstep.historyRegions[rp2key].historyOutputs[ho2key].data
lpf_hist = RIKSstep.historyRegions[asskey].historyOutputs[hoasse].data
maxpos = load_hist.index(max(load_hist,key=lambda x:x[1]))
load = load_hist[maxpos][1]
disp = -disp_hist[maxpos][1]
lpf = lpf_hist[maxpos][1]
odbAccess.closeOdb(myOdb)

# Write the results in the file
out_file.write('\n-RESULTS\n')
out_file.write('Maximum LPF:......................................................... '+str(lpf)+'\n')
out_file.write('Max load:............................................................ '+str(load / 1000)+' [kN]\n')
out_file.write('Displacement at maximum load:........................................ '+str(disp)+' [mm]\n')

# Close the output file
out_file.close()

# Save the cae model
mdb.saveAs(pathName = os.getcwd()+'/'+IDstring+'.cae')

# Return to parent directory
os.chdir('..')

# Write a file where all the results of a batch run are to be gathered
batch_out_file = open('./batch.dat', 'a')
batch_out_file.write("%7.1f %7.1f %3s %4.2f %5.2f %7.3E %6.3f %7.3E %7.3E %7.3E %7.3E %7.3E"%(l_leg, t_shell, int(p_class), lambda_flex, lace_over_leg, l_tot, lpf, load, N_cr_max, N_cr_min, N_cr_tor, N_cr)+'\n')
batch_out_file.close()