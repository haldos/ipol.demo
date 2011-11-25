"""
online image processing
"""

from lib import base_app, build, http, image, config
from lib.misc import app_expose, ctime
from lib.base_app import init_app
import cherrypy
from cherrypy import TimeoutError
import os.path
import shutil

class app(base_app):
	""" template demo app """

	title = "Combined Local Global Optical Flow"
	is_test = True		 # switch to False for deployment
	is_listed = True

	def __init__(self):
		"""
		app setup
		"""
		# setup the parent class
		base_dir = os.path.dirname(os.path.abspath(__file__))
		base_app.__init__(self, base_dir)

		# select the base_app steps to expose
		# index() is generic
		#app_expose(base_app.index)
		#app_expose(base_app.input_select)
		#app_expose(base_app.input_upload)
		# params() is modified from the template
		#app_expose(base_app.params)
		# run() and result() must be defined here

	def build(self):
		"""
		program build/update
		"""
		## store common file path in variables
		tgz_file = self.dl_dir + "clg-87dfe8b.tar.gz"
		prog_file = self.bin_dir + "test_clg_flow"
		log_file = self.base_dir + "build_clg.log"
		## get the latest source archive
		build.download("http://dev.ipol.im/git/?p=haldos/clg.git;a=snapshot;h=87dfe8b1454e0602df8c507a8e7417c3b0a94dcb;sf=tgz", tgz_file)
		## test if the dest file is missing, or too old
		if (os.path.isfile(prog_file)
			 and ctime(tgz_file) < ctime(prog_file)):
			 cherrypy.log("not rebuild needed",
						  context='BUILD', traceback=False)
		else:
			 # extract the archive
			 build.extract(tgz_file, self.src_dir)
			 # build the program
			 build.run("cd %s && cmake . && make" % (self.src_dir + "clg-87dfe8b"), 
						stdout=log_file)
			 print("make complete")
			 # save into bin dir
			 #if os.path.isdir(self.bin_dir):
			 #	 shutil.rmtree(self.bin_dir)
			 #os.mkdir(self.bin_dir)
			 shutil.copy(self.src_dir + os.path.join("clg-87dfe8b", "bin") + "/test_clg_flow", prog_file)
			 print("%s copied to bin directory" % prog_file)
			 # cleanup the source dir
			 shutil.rmtree(self.src_dir)
			 # cleanup the dl dir
			 shutil.rmtree(self.dl_dir)
		return


	@cherrypy.expose
	@init_app
	def wait(self, **kwargs):
		"""
		params handling and run redirection
		"""
		print("ENTER wait")
		print("kwargs = " + str(kwargs))
		self.cfg['param']['alpha'] = kwargs['alpha']
		self.cfg['param']['rho'] = kwargs['rho']
		self.cfg['param']['sigma'] = kwargs['sigma']
		self.cfg['param']['numit'] = kwargs['numit']
		self.cfg.save()
		http.refresh(self.base_url + 'run?key=%s' % self.key)
		return self.tmpl_out("wait.html")

	@cherrypy.expose
	@init_app
	def run(self, **kwargs):
		"""
		algo execution
		"""
		print("ENTER run")
		print("kwargs = " + str(kwargs))
		alpha = self.cfg['param']['alpha']
		rho = self.cfg['param']['rho']
		sigma = self.cfg['param']['sigma']
		numit = self.cfg['param']['numit']
		## run the algorithm
		try:
			 self.run_algo(alpha,rho,sigma,numit)
		except TimeoutError:
			 return self.error(errcode='timeout') 
		except RuntimeError:
			 return self.error(errcode='runtime')
		http.redir_303(self.base_url + 'result?key=%s' % self.key)

		## archive
		##if self.cfg['meta']['original']:
		ar = self.make_archive()
		ar.add_file("a.png", info="input 1")
		ar.add_file("b.png", info="input 2")
		ar.add_file("t.png", info="true flow")
		ar.add_file("stuff_clg.png", info="CLG calculated flow")
		ar.add_info({"alpha": alpha})
		ar.add_info({"rho": rho})
		ar.add_info({"sigma": sigma})
		ar.add_info({"numit": numit})
		ar.save()

		return self.tmpl_out("run.html")


	def run_algo(self, alpha, rho, sigma, numit):
		"""
		the core algo runner
		could also be called by a batch processor
		this one needs no parameter
		"""
		print('ENTERING run_algo')
		print('alpha = ' + str(alpha))
		print('rho = ' + str(rho))
		print('sigma = ' + str(sigma))
		print('numit = ' + str(numit))
		p = self.run_proc(['clgstuff.sh', str(alpha), str(rho), str(sigma), str(numit)])
		self.wait_proc(p, timeout=self.timeout)
		return

	@cherrypy.expose
	def input_select(self, **kwargs):
		"""
		use the selected available input images
		"""
		print("ENTERING input_select")
		self.new_key()
		self.init_cfg()
		print("key = " + self.key)
		# kwargs contains input_id.x and input_id.y
		input_id = kwargs.keys()[0].split('.')[0]
		assert input_id == kwargs.keys()[1].split('.')[0]
		# get the images
		input_dict = config.file_dict(self.input_dir)
		idir = self.input_dir + input_dict[input_id]['subdir'] + '/'
		fnames = "a.png b.png t.png t.tiff".split()
		for f in fnames:
				shutil.copy(idir + f, self.work_dir + f)
		#fnames = input_dict[input_id]['files'].split()
		#for i in range(len(fnames)):
		#	 shutil.copy(self.input_dir + fnames[i],
		#				 self.work_dir + 'input_%i' % i)
		#msg = self.process_input()
		self.log("input selected : %s" % input_id)
		self.cfg['meta']['original'] = False
		self.cfg.save()
		# jump to the params page
		return self.params(msg="(no message)", key=self.key)

	@cherrypy.expose
	def index(self):
		"""
		demo presentation and input menu
		"""
		print("ENTERING index")
		# read the input index as a dict
		inputd = config.file_dict(self.input_dir)
		print(inputd)
		for (input_id, input_info) in inputd.items():
			fname = ["a.png", "b.png", "t.png"]
			inputd[input_id]['baseinput'] = self.input_url + input_info['subdir']+'/'
			inputd[input_id]['url'] = [self.input_url
				  	        + input_info['subdir']
					+ '/' + os.path.basename(f)
			      for f in fname]
			inputd[input_id]['tn_url'] = [self.input_url
						 + input_info['subdir']
						 + '/'
						 + os.path.basename(f)
				 for f in fname]

		return self.tmpl_out("input.html", inputd=inputd)

	@cherrypy.expose
	@init_app
	def result(self):
		"""
		display the algo results
		"""
		return self.tmpl_out("result.html")

# vim: set ts=8 noexpandtab:
