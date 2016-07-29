'''this tool is to calculate the total fluescence yield measured for solid and liquids
the data used is either generated from the installed xraylib http://ftp.esrf.eu/pub/scisoft/xraylib/readme.html
or taken from a database that was previously generated from this source
all set values are handled in three panda Dataframes and explained in them
The frames are stored in a folder that can be either specificly given or the standard folder is used
I can highly recomment to use the xraylib life since additional functions like the compound parser will be available
see the import values function for more details
ALL VALUES USED FOR CALCULATION ARE IN SI!!!!
developed by 
Jens Uhlig 2013
The compund library was snatched from Bruce ravels Hephaestus

FIX SOLVENT CONCENTRATION
FIX GENERAL ENERGY
'''
from __future__ import division
import os,sys,numpy,pylab,re,scipy.constants
import pandas as pd
from pandas import Series,DataFrame
from numpy import exp,log,log10,sin,cos,pi,arange,logspace,absolute
from scipy.interpolate import UnivariateSpline
from matplotlib.ticker import FuncFormatter
import fancy_parser
pylab.ion(); pylab.show()
if 1:#check if there is an load xraylib in the path and attempt to load the x-raylib
	if any('xraylib' in s for s in sys.path):# check if the xraylib is installed and in the python path
		try:
			import _xraylib as xl
			found_xraylib=True
		except:
			raise ImportError('Importing the library went wrong check the installation')
	else:
		found_xraylib=False

if 1: #names of shells and lines and atoms behind if to use code folding Libswitch is here a global that pretty much sets the default for using xraylib or local database
	Atom_Names=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu']
	file_path=os.path.dirname(os.path.realpath(__file__))
	database_path=file_path+os.sep+'databases'
	Libswitch=found_xraylib
	Libswitch=True

def read_lib_entries(from_xraylib=False,regen_Database_from_xraylib=False):
	'''the xraylib has an interesting way to parse the name to the integer that is needed to retrieve the values from the library to avoid recoding the entire project i wrote this little namine parser. It retruns first a Pandas.Series with the shell name as index and as a secon a Pandas.Dataframe with a double index, the first index is from which shell the second index is to which shell the transition goes. The third is a series with pseudo names like Kalpha.'''		
	if from_xraylib:
		names_list=dir(xl)
		lines=[re.split('_',s)[0] for s in names_list if 'LINE' in re.split('_',s)[-1]]
		shells=[re.split('_',s)[0] for s in names_list if 'SHELL' in re.split('_',s)[-1]]
		shells.sort()
		shell_int=[eval('xl.%s_SHELL'%st) for st in shells]
		shells=[s[:2] for s in shells]
		lines_sp=[re.findall(r'([A-Z][.\0-9]*)', l) for l in lines if not ('A' in l or 'B' in l or 'G' in l or 'P' in l or 'O' in l or 'Q' in l or len(re.findall(r'([A-Z][.\0-9]*)', l)[1])<2)]
		lines_int=[eval('xl.%s%s_LINE'%(st[0],st[1])) for st in lines_sp]
		pseudo=[l for l in lines if (('A' in l or 'B' in l) and len(l)>2)]
		pseudo_int=[]
		pseudo_name=[]
		name_dic={'KA1':'K-L3','KA2':'K-L2','KA3':'K-L1','KB1':'K-M3','KB2':'K-N2','KB3':'K-M2','KB4':'K-N4','KB5':'K-M4'}
		for ps in pseudo:
			if ps in name_dic.keys():
				pseudo_name.append(ps)
				pseudo_int.append(name_dic[ps])
		renamed=[]
		for st in pseudo_name:
			b=list(st)
			try:
				b[b.index('A')]='alpha'
			except:
				pass
			try:
				b[b.index('B')]='beta'
			except:
				pass
			renamed.append(''.join(b))
		shell_series=Series(shell_int,index=shells)
		shell_series.index.name='Shell'
		lines_series=Series(lines_int,index=pd.MultiIndex.from_tuples(lines_sp))
		lines_series.index.name=['end_shell','excited_shell']
		pseudo_series=Series(pseudo_int,index=renamed)
		pseudo_series.index.name='pseudoname'
		if regen_Database_from_xraylib:
			lines_series.to_csv(database_path+os.sep+'Line_names.csv')
			shell_series.to_csv(database_path+os.sep+'Shell_names.csv')
			pseudo_series.to_csv(database_path+os.sep+'Pseudo_names.csv')
	else:
		lines_series=Series.from_csv(database_path+os.sep+'Line_names.csv',index_col=[0,1])
		shell_series=Series.from_csv(database_path+os.sep+'Shell_names.csv',index_col=0)
		pseudo_series=Series.from_csv(database_path+os.sep+'Pseudo_names.csv',index_col=0)
		pseudo_series.index.name='pseudoname'
		lines_series.index.name=['end_shell','excited_shell']
		shell_series.index.name='Shell'
	return shell_series,lines_series,pseudo_series
def get_all_shell_names(from_xraylib=False):
	shell_series,_,_=read_lib_entries(from_xraylib=from_xraylib)
	return [str(s) for s in list(shell_series.index)]
def get_all_line_names(from_xraylib=False):
	'''list with all lines'''
	_,line_series,_=read_lib_entries(from_xraylib=from_xraylib)
	return ['%s-%s'%(str(s[0]),str(s[1])) for s in list(line_series.index)]
def get_all_pseudo_names(from_xraylib=False):
	_,_,pseudo_series=read_lib_entries(from_xraylib=from_xraylib)
	return [str(s) for s in list(pseudo_series.index)]
def set_up_data_structure(input=None):
	if input is None:
		zahlen,switches,text=read_standard()
		print 'used saved values'
	else:
		zahlen,switches,text=input
	try:
		za_old=zahlen.copy()
		sw_old=switches.copy()
		st_old=text.copy()
		old_exist=True
	except:
		old_exist=False
		pass
	data_zahlen={	'Field':			['Thickness','Theta','Phi','BeamFlux','Crystals','BeamSize','DetecReflec','DetecLine','DetecRadius','DetecDis','DetecAbsorb','EnIn','EnOut','SamMol','SamDen','Sam2Mol','Sam2Den','SolDen','SolCon','m_fac','SamRel'],
			'description':		['sample thickness','input angle','output angle','incoming photon flux','Number of crystals in the setup','vertical beam size', 'crystal refl. in bragg','spectral fraction of emission line captured','Radius of the Detector','sample Detector distance','Fraction of absorption in space','Excitation energy','emission energy','Sample Molar concentration','Sample Mass density','secon sample Molar con','second sample density','Solvent mass density','Solvent Concentration (rest water)','the m factor is for playing with the self absorption','relation between the samples']}
	data_switch={	'Field':			['xraylib','type','self_absorb','SamRel'],
			'description':		['switch to indicate if xraylib was used','sample type 0-liquid, 1-solid','turns on and of the self absorption correction','relation between the samples']}
	data_text={	'Field':			['exp','Absorber','lines','sample_formula','sample_2_formula','solvent_formula','comment','inner_loop','what_inner','outer_loop','what_outer'],
			'description':		['name of experiment folder','Atom of interest','emission lines, can be list','Formulas, Elements start large, floats allowed','Formulas, Elements start large, floats allowed','same like sample formula','well here we can leave a comment','stores the string of the inner loop','stores what was scanned as inner loop','stores the string of the outer loop','stores what is looped in the outer']}
	zahlen=DataFrame(data_zahlen,columns=['description','value','SIvalue'],index=data_zahlen['Field'])
	zahlen.index.name='Field'
	switches=DataFrame(data_switch,columns=['description','value'],index=data_switch['Field'],dtype=int)
	switches.index.name='Field'
	text=DataFrame(data_text,columns=['description','value'],index=data_text['Field'])
	text.index.name='Field'
	zahlen['SIconvert']=Series([1.e-3,pi/180.,pi/180.,1.e-6,0.001,0.001,0.001,0.001,0.01],index=['Thickness','Theta','Phi','BeamSize','DetecRadius','SamDen','Sam2Den','SolDen','SolCon'])
	zahlen['units']=Series(['mm','degree','degree','micron','mm','m','eV','eV','mMol/l','g/cm3','mMol/l','g/cm3','g/cm3','\%','Photons/s'],index=['Thickness','Theta','Phi','BeamSize','DetecRadius','DetecDis','EnIn','EnOut','SamMol','SamDen','Sam2Mol','Sam2Den','SolDen','SolCon','BeamFlux'])
	if old_exist:
		zahlen['value']=zahlen['value'].combine_first(za_old['value'])
		zahlen['SIconvert']=zahlen['SIconvert'].combine_first(za_old['SIconvert'])
		switches['value']=switches['value'].combine_first(sw_old['value'])
		text['value']=text['value'].combine_first(st_old['value'])
		zahlen.SIconvert=zahlen.SIconvert.fillna(1.)
		zahlen.SIvalue=zahlen.value*zahlen.SIconvert
	switches['value']['xraylib']=Libswitch
	return zahlen,switches,text
def read_standard(exp='standard',print_folder=False):
	'''this program reads the standard file, i use markes in the file to determine where the different sections are stored'''
	path_used = file_path+os.sep+'experiments'+os.sep+exp
	if print_folder:print 'experiment path is: ' + path_used
	if not os.path.exists(path_used):
		print 'path does not exist'
		return False
	zahlen=pd.read_csv(path_used+os.sep+'zahlen.csv',index_col=0)
	zahlen['SIvalue']=zahlen['value']*zahlen.SIconvert
	switches=pd.read_csv(path_used+os.sep+'switches.csv',index_col=0)
	text=pd.read_csv(path_used+os.sep+'text.csv',index_col=0)
	return zahlen,switches,text
def write_standard(Input=None,exp=None):
	'''takes the imput of the three panda dataframes in the order zahlen,switches,text and writes them as comma separated values to disk in the folder given as string after \"exp\", if no experiment is given it dumps it in the folder standard. Everything that was previously stored in this folder will be erased'''
	if input is None:
		raise IOError('you have to give some values to write')
	else:
		zahlen,switches,text=Input
	zahlen['SIvalue']=zahlen['value']*zahlen.SIconvert
	text.value['lines']=Lines_name(text.value['lines'])
	if exp is None:
		if not pd.isnull(text.value['exp']):
			exp=text.value['exp']
		else:
			exp='standard'
	path_used = file_path+os.sep+'experiments'+os.sep+exp
	print 'experiment path is: ' + path_used
	if not os.path.exists(path_used):
		os.makedirs(path_used)
	for a in ['zahlen','switches','text']:
		stri=path_used+os.sep+a+'.csv'
		if os.path.exists(stri):
			os.remove(stri)
		temp=eval(a)
		temp.to_csv(stri)
	return True
def compounds(stringen=None):
	standards=pd.read_csv(database_path+os.sep+'compound.csv',index_col=0)
	if stringen is None:
		listen=list(standards.index)
		for atom in Atom_Names:
			listen.append(atom)
		return listen
	if stringen in standards.index:
		return (standards.ix[stringen]['Formula'],standards.ix[stringen]['Density'])
	elif stringen in standards['Formula'].values:
		indexen=list(standards['Formula'].values).index(stringen)
		return (standards['Formula'].ix[indexen],standards['Density'].ix[indexen])
	else:
		return None
def chemparser(str):
	#return re.findall(r'([A-Z][a-z]*)([.\0-9]*)', str) #that was my parser, but the other one is better because it can parse brackets 
	return fancy_parser.wrapper(str)
def relation_wrapper(compound_1=None,compound_2=None,type_of_relation=None,relation_number=None,Mol_density=None,Mass_density=None):
	'''This little wrapper takes two compounds as strings (sum formulas but can take brakets) it needs a type of relation which is either the string \'wp\'=\'weight percent\' or ppm (parts per million) and a number, as well as either the Mol_density or the Mass_density of compound_2 it returns the Mass_density of compound_1. In the end it calls element masses with the relation applied to the density'''
	if compound_1 is None or compound_2 is None or type_of_relation is None or relation_number is None or (Mol_density is None and Mass_density is None):raise IOError('enter all inforamtion needed (see help)')
	comp_2=elementmasses(compound_2,Mol_density=Mol_density,Mass_density=Mass_density)
	comp2_molar_weight=comp_2.Atomic_weight.mul(comp_2.rel_occ).sum()
	comp_1=elementmasses(compound_1,Mol_density=1.)
	comp1_molar_weight=comp_1.Atomic_weight.mul(comp_1.rel_occ).sum()
	if 'wp' in type_of_relation:#ok we are running weigth percent
		if pd.isnull(Mol_density):#ok we have the massdensity of the second compound
			return Mass_density*relation_number/100
		else:
			return Mol_wrapper(Mol_density*(comp1_molar_weight/comp2_molar_weight)*relation_number/100,comp_1).ele_mass.sum()
	elif 'ppm' in type_of_relation:#we assume we got the parts per million
		if pd.isnull(Mol_density):#ok we have the massdensity of the second compound
			return Mass_density*(comp1_molar_weight/comp2_molar_weight)*relation_number*1e-6
		else:
			return Mol_wrapper(Mol_density*relation_number*1e-6,comp_1).ele_mass.sum()
	else:
		raise IOError('realtion Type not recognized')
def Mol_wrapper(Mol_density,Phys_value):
		Phys_value=Phys_value.copy()
		molecules_per_m3=(Mol_density)*scipy.constants.Avogadro
		mass_constant_in_kg=scipy.constants.atomic_mass
		Phys_value['ele_mass']=Phys_value['rel_occ']*Phys_value['Atomic_weight']*molecules_per_m3*mass_constant_in_kg
		return Phys_value
def elementmasses(str,Mol_density=None,Mass_density=None,SolCon=None):
	'''takes chemical formula and one of the densities in Mol/m3 kg/m3 and returns a pandas DataFrame with element in index and element density kg/m3 ( to kg/m3) as data  THE CONCENTRATION HOOK DOESN'T WORK YET'''
	if pd.isnull(str): return None
	if not pd.isnull(str):comp=compounds(str)
	if not pd.isnull(comp):str,Mass_density=comp
	if not pd.isnull(SolCon):pass
	elements=[x[0] for x in chemparser(str)]
	relative_weights=[x[1] if x[1] is not '' else '1' for x in chemparser(str)]
	relative_weights=[float(s) for s in relative_weights]
	Phys_value=get_phys_values(elements,-1)
	Phys_value['rel_occ']=relative_weights
	if hasattr(Mol_density, '__iter__'):
		Phys_value_list=[]
		for den in Mol_density:
			Phys_value_list.append(Mol_wrapper(den,Phys_value))
		return Phys_value_list
	elif not pd.isnull(Mol_density):		
		Phys_value=Mol_wrapper(Mol_density,Phys_value)
	elif not pd.isnull(Mass_density):
		micro_unit=Mass_density/(sum(Phys_value['rel_occ']*Phys_value['Atomic_weight']))
		Phys_value['ele_mass']=Phys_value['rel_occ']*Phys_value['Atomic_weight']*micro_unit
	else:
		raise IOError('please enter one kind of density or a known compound')
	return Phys_value
def Atoms_number(Atoms,uselist=False):
	'''returns an integer or a list of integers with the atomic numbers of the given atoms'''
	if not hasattr(Atoms, '__iter__'):Atoms=[Atoms]#well lets make it all to a list
	if type(Atoms[0]) is str:#we know now that it is a string
		if re.match(r'[A-Z]',Atoms[0]) is None:#ok the string is a kind of number lets convert this into int
			listen=[int(eval(A)) for A in Atoms]
		else:#ok it is a kind of string so lets thread it like one from the list
			listen=[Atom_Names.index(A)+1 for A in Atoms]
	else:#ok we have already a number
		listen=[int(A) for A in Atoms]
	if not uselist:
		if len(listen) <2:listen=listen[-1]
	return listen
def Atoms_name(Atoms,uselist=False):
	'''returns a string or a list of strings with the name or the Atoms'''
	if not hasattr(Atoms, '__iter__'):Atoms=[Atoms]#well lets make it all to a list
	if type(Atoms[0]) is str:#we know now that it is a string
		if re.match(r'[A-Z]',Atoms[0]) is None:#ok the string is a kind of number les hope the atom number
			listen=[Atom_Names[int(eval(A))-1] for A in Atoms]
		else:# well the atoms are alrady named
			listen=[str(A) for A in Atoms]
	else:#well it is some kind of number
		listen=[Atom_Names[int(A)-1] for A in Atoms]
	if not uselist:
		if len(listen) <2:listen=listen[-1]
	return listen
def Lines_name(Lines,uselist=False):
	'''returns a string or a list of strings with the name or the Lines'''
	_,line_series,_=read_lib_entries(from_xraylib=False)#get a pandas series with all the names as index and the numbers as values
	if not hasattr(Lines, '__iter__'):Lines=[Lines]#well lets make it all to a list
	if type(Lines[0]) is str:#we know now that it is a string
		if re.match(r'[A-Z]',Lines[0]) is None:#ok the string is a kind of number lets hope the Line number
			listen=[]
			for L in Lines:
				a=line_series.index[list(line_series.values).index(int(eval(L)))]
				listen.append('%s-%s'%(str(a[0]),str(a[1])))
		else:# well the Lines are alrady named
			listen=Lines
	else:#ok the stuff are actually numbers
		listen=[]
		for L in Lines:
			a=line_series.index[list(line_series.values).index(L)]
			listen.append('%s-%s'%(str(a[0]),str(a[1])))
	if not uselist:
		if len(listen) <2:listen=listen[-1]
	return listen
def Lines_number(Lines,uselist=False):
	'''returns a string or a list of strings with the numer or the Lines'''
	_,line_series,_=read_lib_entries(from_xraylib=False)#get a pandas series with all the names as index and the numbers as values
	if not hasattr(Lines, '__iter__'):Lines=[Lines]#well lets make it all to a list
	if type(Lines[0]) is str:#oki we have to do with strings
		if re.match(r'[A-Z]',Lines[0]) is None:#ok the string is a kind of number les hope the lines number
			listen=[int(eval(A)) for A in Lines]
		else:# we know there are some letters in it so lets retrieve the index
			listen=[int(line_series[tuple(re.split('-',L))]) for L in Lines]
	else:#well the list is kind of numbers already lets make sure they are all int
		listen=[int(A) for A in Lines]
	if not uselist:
		if len(listen) <2:listen=listen[-1]
	return listen
def Shells_name(Shells,uselist=False):
	'''returns a string or a list of strings with the name or the Shells'''
	shell_series,_,_=read_lib_entries(from_xraylib=False)
	if not hasattr(Shells, '__iter__'):Shells=[Shells]#lets make it all a list
	if type(Shells[0]) is str:#oki we start from a string
		if re.match(r'[A-Z]',Shells[0]) is None:#ok the string is a kind of number les hope the shells number
			listen=[str(shell_series.index[list(shell_series.values).index(int(eval(S)))]) for S in Shells]
		else:# we know there are some letters lets make sure this are normal strings
			listen=[str(S) for S in Shells]
	else:# well the entries are a kind of numbers already
		listen=[str(shell_series.index[list(shell_series.values).index(int(S))]) for S in Shells]
	if not uselist:
		if len(listen) <2:listen=listen[-1]
	return listen
def Shells_number(Shells,uselist=False):
	'''returns a string or a list of strings with the numer or the Lines'''
	shell_series,_,_=read_lib_entries(from_xraylib=False)
	if not hasattr(Shells, '__iter__'):Shells=[Shells]#lets make it all a list
	if type(Shells[0]) is str:#oki we start from a string
		if re.match(r'[A-Z]',Shells[0]) is None:#ok the string is a kind of number les hope the shells number
			listen=[int(eval(S)) for S in Shells]
		else:# we know there are some letters lets find the values in the pandas.series
			listen=[int(shell_series[S[:2]]) for S in Shells]#we only use the first letters, there are a few strange double numbers, deal with this later
	else:# well this are obviously already some numbers, lets make sure these are int
		listen=[int(s) for s in Shells]
	if not uselist:
		if len(listen) <2:listen=listen[-1]
	return listen
def spline_column_wrapper(DatFram,new_x):
	ind=DatFram.index
	if not hasattr(new_x, '__iter__'):new_x=[new_x]
	new_x=['%.3f'%num for num in new_x]
	DatCop=DatFram.copy()
	DatCop=DatCop.reindex(columns=new_x)
	if not hasattr(ind, '__iter__'):ind=[ind]
	for i in ind:
		for j in new_x:
			s = UnivariateSpline(DatFram.columns,DatFram.ix[i], s=0)
			DatCop.ix[i][j]=s(float(j))
	return DatCop
def xraylib_FluorLine_energy_wrapper(Atoms,Line,Energy):
	'''this function is supposed to be called by get_partial_XRF_yield it expects a single atom, a single line, and either a single energy or a energy iterable and returns a pandas series with the energy probed as index'''
	out=[]
	if not hasattr(Energy, '__iter__'):Energy=[Energy]
	for en in Energy:
		if en < xl.LineEnergy(Atoms_number(Atoms),Lines_number(Line)):
			out.append(0.0)
		else:
			out.append(xl.CS_FluorLine_Kissel_Cascade(Atoms_number(Atoms),Lines_number(Line),float(en)/1000.))
	return out
def get_partial_XRF_cross(Atoms=None,Line=None,Energy=None,from_xraylib=Libswitch,regen_Database_from_xraylib=False):
	'''reading and writing the XRF crosssection to file using m2/kg and Kissel Cascade Line index from Xraylib and returning it
	if the switch \"regn_Database_from_xraylib\" is set to True the local Database will be rebuilt
	the switch \"from_xraylib\" triggers if the values are generated from xraylib or are taken fro the local database
	returned is always a pandas DataFrame with the types (names=string or numbers=int) specified in Atoms which can be anything iterable or a single name
	the tricky part is if the database is used and the energy asked for is not in it. Here i call for each element a spline element that fits the existing energy values and returns the spline value between them. This is programmed quite stupid with fixed intervals, so around the edges things get ticky or better the spline smoothes over the edges
	For better results here please use the from_xraylib=True option, they do this more sophisticated'''
	if Atoms is None: raise ValueError('please enter the element for which to retrieve the crossection')
	if Line is None: raise ValueError('please enter the transition line for which to retrieve the crossection')
	if Energy is None: raise ValueError('please enter energy at which to retrieve the crossection')		
	hira1=[]
	hira2=[]
	if not hasattr(Line, '__iter__'):Line=[Line]
	if not hasattr(Atoms, '__iter__'):Atoms=[Atoms]
	if not hasattr(Energy, '__iter__'):Energy=[Energy]
	for atom in Atoms_number(Atoms,uselist=True):
		for line in Lines_number(Line,uselist=True):
			hira1.append(Atoms_name(atom))
			hira2.append(Lines_name(line))
	index_i=pd.MultiIndex.from_arrays([hira1,hira2])
	en_string=[]
	for en in Energy:
		en_string.append('%.3f'%en)
	XRF=DataFrame(columns=en_string,index=index_i)
	XRF.columns.names=['EnIn']
	XRF.index.names=['Elements','Lines']	
	if regen_Database_from_xraylib or from_xraylib:
		for indi in XRF.index:
			if regen_Database_from_xraylib:print indi
			XRF.ix[indi]=xraylib_FluorLine_energy_wrapper(indi[0],indi[1],Energy)  #observe teh libary works in keV
		XRF=XRF*10
		if regen_Database_from_xraylib:
			XRF.to_csv(database_path+os.sep+'XRF.csv')
	else:
		database_XRF=pd.read_csv(database_path+os.sep+'XRF.csv',index_col=[0,1])
		temp=database_XRF.columns
		XRF=spline_column_wrapper(database_XRF.ix[XRF.index],Energy)  #observe teh libary works in keV
		#XRF=spline_column_wrapper(database_XRF.ix[Atoms_name(Atoms)],Energy)
	return XRF
def get_AtomicLevelWidth(Atoms=None,Shells=None,from_xraylib=Libswitch,regen_Database_from_xraylib=False):
	'''reading and writing the Line width from Xraylib and returning it
	if the switch \"regn_Database_from_xraylib\" is set to True the local Database will be rebuilt 
	the switch \"from_xraylib\" triggers if the values are generated from xraylib or are taken fro the local database
	returned is always a pandas DataFrame with the types (names=string or numbers=int) specified in Atoms which can be anything iterable or a single name'''
	if Atoms is None: raise ValueError('please enter the element for which to retrieve the energy')
	if Shells is None: raise ValueError('please enter the Shell for which to retrieve the energy')
	if regen_Database_from_xraylib or from_xraylib:
		hira1=[]
		hira2=[]
		if not hasattr(Shells, '__iter__'):Shells=[Shells]
		if not hasattr(Atoms, '__iter__'):Atoms=[Atoms]
		for atom in Atoms_number(Atoms,uselist=True):
			for shell in Shells_number(Shells,uselist=True):
				hira1.append(Atoms_name(atom))
				hira2.append(Shells_name(shell))
		index_i=pd.MultiIndex.from_arrays([hira1,hira2])
		out=[]
		for ind in index_i:
			out.append(xl.AtomicLevelWidth(Atoms_number(ind[0]),Shells_number(ind[1])))
		AtomicLevelWidth=Series(out,index=index_i)*1000
		AtomicLevelWidth.index.names=['Elements','Shells']	                   
		if regen_Database_from_xraylib:
			AtomicLevelWidth.to_csv(database_path+os.sep+'AtomicLevelWidth.csv')
	else:
		AtomicLevelWidth=pd.read_csv(database_path+os.sep+'AtomicLevelWidth.csv',index_col=[0,1])
		AtomicLevelWidth.index.names=['Elements','Shells']
	return AtomicLevelWidth
def get_LineWidth(Atoms=None,Lines=None,from_xraylib=Libswitch):
	if Atoms is None: raise ValueError('please enter the element for which to retrieve the width')
	if Lines is None: raise ValueError('please enter the Lines for which to retrieve the width')
	#ok need to retrieve level width first split the names and make a list of the levels needed
	shell_list=[]
	transition_list=[]
	for line in Lines_name(Lines,uselist=True):
		help=re.split('-',line)
		transition_list.append(help)
		if help[0] not in shell_list:
			shell_list.append(help[0])
		if help[1] not in shell_list:
			shell_list.append(help[1])
	Levels_width=get_AtomicLevelWidth(Atoms=Atoms,Shells=shell_list)
	Line_energy=get_LineEnergy(Atoms=Atoms,Lines=Lines)
	Line_info={'LineEnergy':Line_energy,'LineWidth':0.0}
	Line_info=DataFrame(Line_info,index=Line_energy.index)

	for atom in Atoms_name(Atoms,uselist=True):	
		for i,line in enumerate(Lines_name(Lines,uselist=True)):
			Line_info['LineWidth'][atom,line]=Levels_width.ix[atom][transition_list[i][0]]+Levels_width.ix[atom][transition_list[i][1]]
	return Line_info
def xraylib_TotCross_energy_wrapper(Atoms,Energy):
	'''this function is supposed to be called by get_total_cross_yield it expects a single atom, a single line, and either a single energy or a energy iterable and returns a pandas series with the energy probed as index'''
	out=[]
	for en in Energy:
		out.append(xl.CS_Total(Atoms_number(Atoms),float(en)/1000.))#the database is in keV
	return out
def get_total_cross(Atoms=None,Energy=None,from_xraylib=Libswitch,regen_Database_from_xraylib=False):
	'''reading and writing the total crosssection to file using m2/kg from Xraylib and returning it
	if the switch \"regn_Database_from_xraylib\" is set to True the local Database will be rebuilt
	the switch \"from_xraylib\" triggers if the values are generated from xraylib or are taken fro the local database
	returned is always a pandas DataFrame with the types (names=string or numbers=int) specified in Atoms which can be anything iterable or a single name
	the tricky part is if the database is used and the energy asked for is not in it. Here i call for each element a spline element that fits the existing energy values and returns the spline value between them. This is programmed quite stupid with fixed intervals, so around the edges things get ticky or better the spline smoothes over the edges
	For better results here please use the from_xraylib=True option, they do this more sophisticated'''
	if Atoms is None: raise ValueError('please enter the element for which to retrieve the crossection')
	if Energy is None: raise ValueError('please enter energy at which to retrieve the crossection')
	if not hasattr(Atoms, '__iter__'):Atoms=[Atoms]
	if not hasattr(Energy, '__iter__'):Energy=[Energy]
	hira1=[]
	for atom in Atoms_number(Atoms,uselist=True):
			hira1.append(Atoms_name(atom))
	en_string=[]
	for en in Energy:
		en_string.append('%.3f'%en)
	TotCross=DataFrame(columns=en_string,index=hira1)
	TotCross.columns.names=['EnIn']
	TotCross.index.names=['Elements']	
	if regen_Database_from_xraylib or from_xraylib:
		for indi in TotCross.index:
			TotCross.ix[indi]=xraylib_TotCross_energy_wrapper(indi,Energy) #observe teh libary works in keV
		TotCross=TotCross*10. #we work SI units
		if regen_Database_from_xraylib:
			TotCross.to_csv(database_path+os.sep+'TotCross.csv')
	else:
		database_TotCross=pd.read_csv(database_path+os.sep+'TotCross.csv',index_col=[0])
		temp=database_TotCross.columns
		TotCross=spline_column_wrapper(database_TotCross.ix[Atoms_name(Atoms)],Energy)#this stuff is fromt he database so this is the right unit
		TotCross.columns.names=['EnIn']
		TotCross.index.names=['Elements']
	return TotCross
def get_absorb(compound=None,Energy=None,density=None,from_xraylib=Libswitch):
	'''get the absorption length in m, enter a Energy in eV and a density in g/cm3, many values are tabulated
	please enter either a compound as string or a valid formula. If you enter a formula you have to provide a density.
	You can check what is in the database by calling the empty function compounds()'''
	if density: density=density*1e3#ok density was given, since i use Kg/m3 i need to change here
	if compound is None:
		print 'please enter a named compound or a formula'
		return False
	if Energy is None:
		print 'please enter a single energy or a list/vector or energies'
		return False
	if not hasattr(Energy, '__iter__'):Energy=[Energy]
	if density: #ok somebody entered a density so we use exactly what is given 
		formula=compound
		density=density
	elif compound in compounds(): #ok we have this in the standards list
		formula,density=compounds(compound)
	else:
		print 'compound name not found in List, please enter either a compound or a valid formula. If you enter a formula you have to provide a density. You can check what is in the database by calling the empty function compounds()'
		return False
		
	df=pd.DataFrame()
	form=chemparser(formula)
	total_mass=0
	for comb in form:
		cross=get_total_cross(Atoms=comb[0],Energy=Energy,from_xraylib=from_xraylib)#atomar crossection
		atomar_mass=get_phys_values(comb[0])['Atomic_weight']*comb[1] #molar mass and how many atoms
		cross=cross*atomar_mass*density
		total_mass+=atomar_mass
		df=	df.append(cross)
	return 100/(df.sum()/total_mass)#I must screw somewhere the numbers up I am a factor 100 wrong but he letters are right. so a cheat here.
def get_EdgeEnergy(Atoms=None,Shells=None,from_xraylib=Libswitch,regen_Database_from_xraylib=False):
	'''reading and writing the Edge Energy from Xraylib and returning it
	if the switch \"regn_Database_from_xraylib\" is set to True the local Database will be rebuilt 
	the switch \"from_xraylib\" triggers if the values are generated from xraylib or are taken fro the local database
	returned is always a pandas DataFrame with the types (names=string or numbers=int) specified in Atoms which can be anything iterable or a single name'''
	if Atoms is None: raise ValueError('please enter the element for which to retrieve the energy')
	if Shells is None: raise ValueError('please enter the Shell for which to retrieve the energy')
	hira1=[]
	hira2=[]
	if not hasattr(Shells, '__iter__'):Shells=[Shells]
	if not hasattr(Atoms, '__iter__'):Atoms=[Atoms]
	for atom in Atoms_number(Atoms,uselist=True):
		for shell in Shells_number(Shells,uselist=True):
			hira1.append(Atoms_name(atom))
			hira2.append(Shells_name(shell))
	index_i=pd.MultiIndex.from_arrays([hira1,hira2])
	out=[]
	if regen_Database_from_xraylib or from_xraylib:
		for ind in index_i:
			out.append(xl.EdgeEnergy(Atoms_number(ind[0]),Shells_number(ind[1])))
		EdgeEnergy=Series(out,index=index_i)*1000.
		EdgeEnergy.index.names=['Elements','Shells']	                   
		if regen_Database_from_xraylib:
			EdgeEnergy.to_csv(database_path+os.sep+'EdgeEnergy.csv')
	else:
		readEdgeEnergy=pd.read_csv(database_path+os.sep+'EdgeEnergy.csv',index_col=[0,1])
		for ind in index_i:
			out.append(readEdgeEnergy.ix[ind])
		EdgeEnergy=Series(out,index=index_i)
		EdgeEnergy.index.names=['Elements','Shells']	
	return EdgeEnergy
def get_LineEnergy(Atoms=None,Lines=None,from_xraylib=Libswitch,regen_Database_from_xraylib=False):
	'''reading and writing the Edge Energy from Xraylib and returning it
	if the switch \"regn_Database_from_xraylib\" is set to True the local Database will be rebuilt 
	the switch \"from_xraylib\" triggers if the values are generated from xraylib or are taken fro the local database
	returned is always a pandas DataFrame with the types (names=string or numbers=int) specified in Atoms which can be anything iterable or a single name'''
	if Atoms is None: raise ValueError('please enter the element for which to retrieve the energy')
	if Lines is None: raise ValueError('please enter the Shell for which to retrieve the energy')
	hira1=[]
	hira2=[]	
	if not hasattr(Lines, '__iter__'):Lines=[Lines]
	if not hasattr(Atoms, '__iter__'):Atoms=[Atoms]
	for atom in Atoms_number(Atoms,uselist=True):
		for line in Lines_number(Lines,uselist=True):
			hira1.append(Atoms_name(atom))
			hira2.append(Lines_name(line))
	index_i=pd.MultiIndex.from_arrays([hira1,hira2])
	out=[]
	if regen_Database_from_xraylib or from_xraylib:
		for ind in index_i:
			out.append(xl.LineEnergy(Atoms_number(ind[0]),Lines_number(ind[1])))
		LineEnergy=Series(out,index=index_i)*1000#we think eV
		LineEnergy.index.names=['Elements','Lines']	                   
		if regen_Database_from_xraylib:
			LineEnergy.to_csv(database_path+os.sep+'LineEnergy.csv')
	else:
		read_LineEnergy=pd.read_csv(database_path+os.sep+'LineEnergy.csv',index_col=[0,1])
		for ind in index_i:
			out.append(read_LineEnergy.ix[ind])
		LineEnergy=Series(out,index=index_i)
		LineEnergy.index.names=['Elements','Lines']
	return LineEnergy
def get_phys_values(Atoms=None,Shells=None,from_xraylib=Libswitch,regen_Database_from_xraylib=False):
	'''getting the atomic data like weight name and shell occupation and returning it
	if the switch \"regn_Database_from_xraylib\" is set to True the local Database will be rebuilt
	the switch \"from_xraylib\" triggers if the values are generated from xraylib or are taken fro the local database
	returned is always a pandas DataFrame with the types (names=string or numbers=int) specified in Atoms which can be anything iterable or a single name'''
	if Atoms is None: raise ValueError('please enter the element for which to retrieve the energy')
	if Shells is None: Shells=Shells_number(get_all_shell_names())
	if regen_Database_from_xraylib or from_xraylib:
		Phys_dat={}
		Phys_dat['Atomic_number']=Atoms_number(Atoms)
		weight=[]
		if not hasattr(Atoms, '__iter__'):
			Atoms=[Atoms]
		elif type(Atoms) is str:
			Atoms=[Atoms]
		for x in Atoms_number(Atoms,uselist=True):
			weight.append(xl.AtomicWeight(x))
		Phys_dat['Atomic_weight']=weight
		Atom=DataFrame(Phys_dat,columns=('Atomic_number','Atomic_weight'),index=Atoms_name(Atoms,uselist=True))
		if Shells is not -1:
			Ele_conf=DataFrame(columns=Shells_name(Shells),index=Atoms_name(Phys_dat['Atomic_number'],uselist=True))
			for x in Atoms_number(Phys_dat['Atomic_number'],uselist=True):
				temp=[]
				for i in Shells_number(Shells,uselist=True):
					temp.append(xl.ElectronConfig(x,i))
				Ele_conf.ix[Atoms_name(x)]=temp
			Atom_phys=pd.merge(Atom,Ele_conf,left_index=True,right_index=True)
		else:
			Atom_phys=Atom.copy()
		Atom_phys.index.name='Elements'
		if regen_Database_from_xraylib:
			Atom_phys.to_csv(database_path+os.sep+'Atom_phys.csv')
		if from_xraylib:
			return Atom_phys.ix[Atoms_name(Atoms)]
	else:
		Atom_phys=pd.read_csv(database_path+os.sep+'Atom_phys.csv',index_col=0)
	return Atom_phys.ix[Atoms_name(Atoms)]	
def rebuild_complete_database(tousandsofenergy=20):
	raw_input("This will take about 15h to rebuild the complete database, press \"enter\" to run or \"ctrl-c\" to abort ")
	#EnIn_range=range(1000,20000,energystep)
	Atoms=Atom_Names
	shell,lines,_=read_lib_entries(from_xraylib=True,regen_Database_from_xraylib=True)
	print 'Lists rebuild'
	Line=['%s-%s'%(a[0],a[1]) for a in lines.index]
	Shells=Shells_number(list(shell.index))
	#Atoms=Atoms_name([25,26])#range(10,28,1)#Atom_Names
	# Line=range(0,3,1)#Line_Names
	# Shells=range(0,3,1)	
	_=get_phys_values(Atoms=Atoms,Shells=Shells,regen_Database_from_xraylib=True)
	print 'Phys_values rebuild'
	temp=get_EdgeEnergy(Atoms=Atoms,Shells=Shells,regen_Database_from_xraylib=True)
	#Lets build a list of energies for which we actually want to make all the tables
	EnIn=[1000]
	while EnIn[-1]<temp.max():
		EnIn.append(EnIn[-1]+EnIn[-1]*tousandsofenergy/1000.)
	around_edge=[-0.1,0,0.1]
	for en in temp.values:
		if en > 1000:
			for d_en in around_edge:
				EnIn.append(en+d_en)
	EnIn.sort()
	EnIn_range=[1000]
	for s in EnIn:
		if s > EnIn_range[-1]+0.2:
			EnIn_range.append(s)
	print 'Edge energies rebuild'	
	_=get_AtomicLevelWidth(Atoms=Atoms,Shells=Shells,regen_Database_from_xraylib=True)
	print 'Atomic Level width rebuild'
	_=get_LineEnergy(Atoms=Atoms,Lines=Line,regen_Database_from_xraylib=True)
	print 'Line Energies rebuild width rebuild'
	_=get_total_cross(Atoms=Atoms,Energy=EnIn_range,regen_Database_from_xraylib=True)
	print 'total crossections rebuild'
	_=get_partial_XRF_cross(Atoms=Atoms,Line=Line,Energy=EnIn_range,regen_Database_from_xraylib=True)
	print 'Partial crossections rebuild'
def emission_detection(Input=None,EnIn=None,EnOut=None,SamMol=None,Thickness=None,m_fac=None):
	if Input is None:
		(zahlen,switches,text)=read_standard()
	else:
		zahlen,switches,text=Input
		zahlen.SIvalue=zahlen.value*zahlen.SIconvert
	za=zahlen['SIvalue'];sw=switches['value'];te=text['value']
	if not sw['self_absorb']:
		za['m_fac']=1.0
		m_fac=1
	varied={}
	if EnIn is None:EnIn=za['EnIn']#this will be handled inside the subfunctions and retur different columns for differnet energies
	else:varied['EnIn']=EnIn
	if EnOut is None:
		EnOut=za['EnOut']#these will be handled below and generate differnet columns
	else:# ok we have a input is it numbers or strings?
		if not hasattr(EnOut, '__iter__'):
			temp=[EnOut]
		if type(EnOut[0]) is str:#somebody gives a list of strings with numbers or names lets attach the energies only
			EnOut_new=[]
			for i in EnOut:
				EnOut_new.append(get_LineEnergy(Atoms=te['Absorber'],Lines=Lines_number(i),from_xraylib=Libswitch))
		varied['EnOut']=EnOut
	if SamMol is None:SamMol=za['SamMol']#this will makea list for the differnet values
	else:varied['SamMol']=SamMol
	if Thickness is None:Thickness=za['Thickness'] #these will make a list for the differnet values
	else:varied['Thickness']=Thickness
	if m_fac is None:m_fac=za['m_fac'] #these will make a list for the differnet values
	else:varied['m_fac']=m_fac
	if len(varied) > 1:raise IOError('Not more then One value can be varied at any time in this function')
	#NEED TO FIX SOLVENT CONCENTRATION HOOK
	if 'None' in te['solvent_formula']:
		solvent=None
	elif pd.isnull(te['solvent_formula']):
		solvent=None
	elif len(te['solvent_formula']) <1:
		solvent=None		
	else:
		solvent=elementmasses(te['solvent_formula'],Mass_density=za['SolDen'],SolCon=za['SolCon'])
	sample=elementmasses(te['sample_formula'],Mol_density=SamMol,Mass_density=za['SamDen'])
	if ('None' in te['sample_2_formula']) or pd.isnull(te['sample_2_formula']):
		sample2=None
	else:
		sample2=elementmasses(te['sample_2_formula'],Mol_density=za['Sam2Mol'],Mass_density=za['Sam2Den'])
	if not te['Absorber'] in sample.index: raise IOError('Absorber not included in sample')
	if not hasattr(SamMol, '__iter__'):#ok we have only one density for the sample
		if sample2 is None:
			if solvent is not None:
				total=sample.ele_mass.add(solvent.ele_mass, fill_value=0)
			else:
				total=sample.ele_mass
		else: #this catches if we actually have two materials (should solve the solid/liquid case
			if solvent is not None:
				total=sample.ele_mass.add(solvent.ele_mass,fill_value=0)
			else:
				total=sample.ele_mass
			total=total.add(sample2.ele_mass,fill_value=0)
	else:# we have a list of densities for the sample
		total=[]
		for i in sample: #loop over the different concentrations
			if sample2 is None:
				total.append(i.ele_mass.add(solvent.ele_mass, fill_value=0))
			else:	#this catches if we actually have two materials (should solve the solid/liquid case
				total_help=i.ele_mass.add(solvent.ele_mass,fill_value=0)
				total_help=total.add(sample2.ele_mass,fill_value=0)
				total.append(total_help)
	#get the absorption cross sections for the incoming energies
	if 1: #total energy munging put in "if" to use code folding
		all_E=[]
		if not hasattr(EnIn, '__iter__'):
			run_EnIn=[EnIn]
		else:
			run_EnIn=EnIn
		for en in run_EnIn:#
			all_E.append(float(en))
		if not 'EnOut' in varied.keys():#check if outgoing energy is given if not generate from element and lines
			if pd.isnull(te['lines']) and pd.isnull(te['Absorber']) and pd.isnull(za['EnOut']):raise IOError('either and energy or emission lines have to be given')
			if not pd.isnull(te['lines']) and not pd.isnull(te['Absorber']):
				Out=get_LineEnergy(Atoms=te['Absorber'],Lines=te['lines'],from_xraylib=Libswitch)
				EnOut=Out
			else:
				EnOut=za['EnOut']
		if not hasattr(EnOut, '__iter__'):
			run_EnOut=[EnOut]
		else:
			run_EnOut=EnOut
		for en in run_EnOut:#add the emission energies to the list to make onle one inquery
			all_E.append(float(en))
	#print all_E
	if not hasattr(SamMol, '__iter__'):#get the absorption crossection for ll elements
		mu=get_total_cross(total.index,all_E)
	else:
		mu=get_total_cross(total[0].index,all_E)
	##############emission#####################
	#print mu
	captured_solid_angle=(za['Crystals']*pi*za['DetecRadius']**2)/(4*pi*za['DetecDis']**2)
	captured_intensity=za['BeamFlux']*za['DetecReflec']*za['DetecAbsorb']*captured_solid_angle
	if hasattr(Thickness, '__iter__'):# ok we vary the samplethickness means the concentration and the EnIn is constant return one chi per
		I_out=[]
		for d in Thickness:
			mu['%.3f'%EnIn][te['Absorber']]=mu['%.3f'%EnIn][te['Absorber']]*m_fac
			chi_rho_d=mu['%.3f'%EnIn].mul(total).sum()*d/sin(za['Theta'])+mu['%.3f'%EnOut].mul(total).sum()*d/sin(za['Phi'])
			Q=get_partial_XRF_cross(Atoms=te['Absorber'],Line=te['lines'],Energy=EnIn)
			q=Q.ix[Q.index[0]]['%.3f'%EnIn]*m_fac
			corr=(1-exp(-chi_rho_d))/chi_rho_d
			I_d=captured_intensity*q*d*total[te['Absorber']]*corr
			I_out.append(I_d)
	elif hasattr(SamMol, '__iter__'):#ok we have only one thickness but more then one concentration? retrun one chi per
		I_out=[]
		for tot in total:
			mu['%.3f'%EnIn][te['Absorber']]=mu['%.3f'%EnIn][te['Absorber']]*m_fac
			chi_rho_d=mu['%.3f'%EnIn].mul(tot).sum()*Thickness/sin(za['Theta'])+mu['%.3f'%EnOut].mul(tot).sum()*Thickness/sin(za['Phi'])
			Q=get_partial_XRF_cross(Atoms=te['Absorber'],Line=te['lines'],Energy=EnIn)
			corr=(1-exp(-chi_rho_d))/chi_rho_d
			q=Q.ix[Q.index[0]]['%.3f'%EnIn]*m_fac
			I_d=captured_intensity*q*Thickness*tot[te['Absorber']]*corr
			I_out.append(I_d)
	elif hasattr(EnIn, '__iter__'):#ok we have only one thickness and one concentration but more then one Excitation energy? retrun one chi per
		I_out=[]
		for en in EnIn:
			mu['%.3f'%en][te['Absorber']]=mu['%.3f'%en][te['Absorber']]*m_fac
			chi_rho_d=mu['%.3f'%en].mul(total).sum()*Thickness/sin(za['Theta'])+mu['%.3f'%EnOut].mul(total).sum()*Thickness/sin(za['Phi'])
			Q=get_partial_XRF_cross(Atoms=te['Absorber'],Line=te['lines'],Energy=en)
			q=Q.ix[Q.index[0]]['%.3f'%en]*m_fac
			corr=(1-exp(-chi_rho_d))/chi_rho_d
			I_d=captured_intensity*q*Thickness*total[te['Absorber']]*corr
			I_out.append(I_d)
	elif hasattr(m_fac, '__iter__'):#ok we have only one thickness, concentration, Excitation energy but more then one m? retrun one chi per
		I_out=[]
		for i in m_fac:
			mu['%.3f'%EnIn][te['Absorber']]=mu['%.3f'%EnIn][te['Absorber']]*i
			chi_rho_d=mu['%.3f'%EnIn].mul(total).sum()*Thickness/sin(za['Theta'])+mu['%.3f'%EnOut].mul(total).sum()*Thickness/sin(za['Phi'])
			Q=get_partial_XRF_cross(Atoms=te['Absorber'],Line=te['lines'],Energy=EnIn)
			q=Q.ix[Q.index[0]]['%.3f'%EnIn]*i
			corr=(1-exp(-chi_rho_d))/chi_rho_d
			I_d=captured_intensity*q*Thickness*total[te['Absorber']]*corr
			I_out.append(I_d)
	else:# one thickness, one concentration, one EnIn,
		mu['%.3f'%EnIn][te['Absorber']]=mu['%.3f'%EnIn][te['Absorber']]*m_fac
		chi_rho_d=mu['%.3f'%EnIn].mul(total).sum()*Thickness/sin(za['Theta'])+mu['%.3f'%EnOut].mul(total).sum()*Thickness/sin(za['Phi'])
		Q=get_partial_XRF_cross(Atoms=te['Absorber'],Line=te['lines'],Energy=EnIn)
		corr=(1-exp(-chi_rho_d))/chi_rho_d
		q=Q.ix[Q.index[0]]['%.3f'%EnIn]*m_fac#assume one EnOut easy to change from here
		I_out=captured_intensity*q*Thickness*total[te['Absorber']]*corr			
	if len(varied) < 1:
		return I_out
	else:
		varied['I_f']=I_out
		return varied
def vary_something(Input=None,what=None,listen=None,is_inner_loop=None):
	'''this is a tiny wrapper fuction to use the varied input
	with what you choose the type, recognized is \"mol\" \"Thickness\" \"EnIn\" and \"m_fac\"as the what variable
	Mol in mMol/l D in mm and EnIn in eV the \"listen\"variable should be a python touple in the manner: (from,to,step)'''
	if listen is None:
		if what is not None:raise IOError('please enter range')
	if what is None and	is_inner_loop is None: raise IOError('please enter something to vary')
	if Input is None:
		zahlen,switches,text=read_standard()
	else:
		zahlen,switches,text=Input
	if 'mol' in what:
		listen=listen*zahlen.SIconvert.ix['SamMol']
		varied=emission_detection(Input=(zahlen,switches,text),SamMol=listen)
	elif 'Thickness' in what:
		listen=listen*zahlen.SIconvert.ix['Thickness']
		varied=emission_detection(Input=(zahlen,switches,text),Thickness=listen)
	elif 'EnIn' in what:
		varied=emission_detection(Input=(zahlen,switches,text),EnIn=listen)
	elif 'm_fac' in what:
		varied=emission_detection(Input=(zahlen,switches,text),m_fac=listen)
	else:
		if Input is not None:# ok we passed a changing input to this and loop outside this funtion
			varied={}
			varied['I_f']=emission_detection(Input=(zahlen,switches,text))
			varied['Nichts']=1
			return varied
		else: raise IOError('something needs to change')
	if is_inner_loop is not None:
		return varied
	for key in varied.keys():
		if 'I_f' in key:
			pass
		else:
			fig = pylab.figure()
			ax = fig.add_subplot(111)
			if 'm' in key:
				plotenlargex=0.05
				plotenlargey=0.05
				ax.plot(varied[key],varied['I_f'],'o')
				middley=emission_detection()
				x=[min(varied[key])-plotenlargex,max(varied[key])+plotenlargex]
				y=[middley+middley*(min(x)-1),middley+middley*(max(x)-1)]
				ax.plot(x,y,'-.')
				pylab.xlim(min(x),max(x))
				pylab.ylim(min(y),max(y))
			else:
				ax.plot(varied[key]/zahlen.SIconvert.ix[key],varied['I_f'])
			ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%g')%(x)))
			ax.set_ylabel('Emission intensity in Photons/s')
			ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%g'%x))
			ax.set_xlabel('%s in %s' %(key,zahlen.units['%s'%key]))
			pylab.title('Emission Intensity vs %s' %key)
	return varied
def vary_a_second(Input=None,what_out=None,what_in=None,listen_in=None,listen_out=None):
	'''this is a tiny wrapper fuction to use the varied input the outer is the one running here and loops over the inner change, the outer change here changes can be done in "zahlen" variable that is handed to vary_something'''
	if listen_in is None and listen_out is None:raise IOError('please enter ranges')
	if what_out is None:raise IOError('please enter what to change')
	if Input is None:
		Input=read_standard()
	zahlen,switches,text=Input
	fig = pylab.figure()
	ax = fig.add_subplot(111)
	pylab.hold('on')
	x_help=0
	ylim=[]
	templist=[]
	for i in listen_out:
		x_help+=1
		if what_out in zahlen.index:
			i_inp=i*zahlen.SIconvert.ix['%s'%what_out]			
			zahlen.value['%s'%what_out]=i_inp
		elif what_out in switches.index:
			switches.value['%s'%what_out]=i
		elif what_out in text.index:
			text.value['%s'%what_out]=i
		else:
			raise IOError('whatout is not a key')
		temp=vary_something(Input=(zahlen,switches,text),what=what_in,listen=listen_in,is_inner_loop=True)
		templist.append(temp)
		for key in temp.keys():
			try:
				ylim.append(temp['I_f'].max())
				ylim.append(temp['I_f'].min())
			except:
				ylim.append(max(temp['I_f']))
				ylim.append(min(temp['I_f']))
				
			if key is 'I_f':
				continue
			elif key is 'm_fac':
				plotenlargex=0.05
				plotenlargey=0.05
				ax.plot(temp[key],temp['I_f'],'o',label='%f'%i)
				middley=emission_detection()
				x=[min(temp[key])-plotenlargex,max(temp[key])+plotenlargex]
				y=[middley+middley*(min(x)-1),middley+middley*(max(x)-1)]
				ax.plot(x,y,'-.',label='%f'%i)
				pylab.xlim(min(x),max(x))
				pylab.ylim(min(y),max(y))
			else:
				if key in zahlen.index:
					ax.plot(temp[key]/zahlen.SIconvert.ix[key],temp['I_f'],'-',label='%s'%i)
				elif 'Nichts' in key:
					ax.plot([x_help],temp['I_f'],'-',label='%s'%i)
				else:
					raise IOError('no key found')
	ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%g'%x))
	if key in zahlen.units:
		ax.set_xlabel('%s in %s' %(key,zahlen.units.ix['%s'%key]))
	ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%g')%(x)))
	ax.set_ylabel('Emission intensity in Photons/s')
	pylab.title('Emission Intensity vs %s' %key)
	ax.legend(loc="best",title='we changed %s'%what_out)
	if min(ylim)>0:
		blub=min(ylim)
	else:
		blub=-max(ylim)/50
	pylab.ylim(blub,max(ylim))
	pylab.hold('off')
	# if abs(max(ylim)-min(ylim)) > 1e5:
		# pylab.yscale('log')
	#ax.legend(loc="center left", bbox_to_anchor=[0.5, 0.5],ncol=2, shadow=True, title="Legend")
	return templist
def run_loop_from_save(Input=None,output=True):
	'''This one is just a little wrapper to make the loop, Output=True returns what was generated to plot'''
	if Input is None:
		Input=read_standard()
	zahlen,switches,text=Input
	
	if pd.isnull(text.value['what_outer']) or 'None' in text.value['what_outer']:#ok we running only inner loop
		varied=vary_something(Input,what=text.value['what_inner'],listen=eval(text.value['inner_loop']),is_inner_loop=None)
	else:
		varied=vary_a_second(Input=Input,what_out=text.value['what_outer'],what_in=text.value['what_inner'],listen_in=eval(text.value['inner_loop']),listen_out=eval(text.value['outer_loop']))
	if output:
		return varied
		
		
if __name__=="__main__":
	if 0:#rebuild complete databsse
		rebuild_complete_database()
	elif 0:#emission detection
		emission_detection()
	elif 0:#vary stuff
		#vary_a_second(what_out='Thickness',listen_out=arange(0.1,1.,0.1))
		vary_a_second(what_in='EnIn',what_out='SamMol',listen_in=arange(7111.95,7112.05,0.001),listen_out=logspace(log10(10),log10(1000),10))
		#vary_something(what='m',listen=arange(0.9,1.1,0.1))
	elif 0:
		zahlen,switches,text=read_standard()
		relation_wrapper('Fe','C','wp',1,Mass_density=3*zahlen.SIconvert['SamDen'])
		print relation_wrapper('Fe','C','ppm',1e6,Mass_density=1*zahlen.SIconvert['SamDen'])/zahlen.SIconvert['SamDen']
	else:
		pass
