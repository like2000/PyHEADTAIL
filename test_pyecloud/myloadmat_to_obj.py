import scipy.io as sio
import numpy as np

def myloadmat(filename, squeeze = True):
	dict_var=sio.loadmat(filename)
	if squeeze:
		for kk in dict_var.keys():
			try:
				dict_var[kk]=np.squeeze(dict_var[kk])
			except:
				pass
	return dict_var
	
class obj_from_dict:
	def __init__(self, dictto):
		for kk in dictto.keys():
			exec 'self.'+kk +'= dictto[kk]'
			
			
def obj_to_dict(obj):
	dict_out={}
	members = dir(obj)
	for member in members:
		exec "dict_out['%s'] = obj.%s"%(member,member)
	return dict_out
			
def myloadmat_to_obj(filename, squeeze = True):
	return  obj_from_dict(myloadmat(filename, squeeze=squeeze))

