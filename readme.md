# HOW-TO use the framework

First you need to define the material object. 

For EELS use `VarName = Material('MaterialName',InitialEnergy);`  

	Mat = Material('Au',2000);
	
Fox XPS use `VarName = XPSMaterial('MaterialName',ShellNumber,'AnodeName');`. `AnodeName` is optional, Al by default.

	Mat = XPSMaterial('Au',1,'Al');
  
  
In the case of Au: 1 shell is 4s1/2, 2 - 4p1/2, 3 - 4p3/2 etc.
Aluminium anode is the default value, so you can write just like this:

`Mat = XPSMaterial('Au',1);`

It is also possible to write this: `VarName = XPSMaterial('MaterialName',ShellNumber,AnodeEnergy);` for example: 

`Mat = XPSMaterial('Au',1,1486.6);`

## Calculations:

There are three possible methods to calculate energy spectra:
* Numerical Solution NS
* Small Angle Approximation SA
* Straight Line Approximation SLA

First it is necessary to set a layer, whether it's finite or infinite:
`Layer(MaterialVarName,thickness in nm);`
* Finite layer: `Layer(Mat,5)`;
* Infinite layer: `Layer(Mat);` or `Layer(Mat,inf);`

#### Reflection (REELS) calculaion
```
VarName = NSReflection(Layer(Mat,5)); %REELS calculation of layer with thickness of 5 nm
	  SAReflection(Layer(Mat));
	  SLAReflection(Layer(Mat));
```
In this class you may set the next:
```
VarName.theta0 = 0; %The incident polar angle with respect to the surface normal
VarName.phi = 0; %The azimuthal angle, 0 in default
VarName.N_in = 15; %The number of inelastic scatterings
VarName.Calculate; %To calculate general charasteristics
```

Then you can calculate the next:
```
VarName.CalculateAngularDistribution;
VarName.CalculateEnergyDistribution(detecting polar angle);
VarName.CalculateEnergyDistribution(45);
VarName.CalculateIneasticScatteringDistribution(detecting polar angle);
```

#### Transmission (TEELS) calculation
```
VarName = NSTransmission(Layer(Mat,5)); %in this case a layer always should be finite (because we calculate transmission through the sample)
```
All other things are the same (like in the reflection case).

#### XPS calculation
```
VarName = NSXPS(Layer(Mat));
````
   All other things are the same (like in the reflection case).

To plot:
```
VarName.plotAngularDistribution;
VarName.plotEnergyDistribution(45);
VarName.plotInelasticScatteringDistribution(45);
```


## To calculate a multilayer system (for example surface and bulk):
### EELS

First it is necessary to set number of layers and their thicknesses.
```
Layers = [Layer(Mat,3) Layer(Mat)]; %if we have just surface and bulk for one material
% or
Layers = [Layer(Mat1,3) Layer(Mat2)]; %if we have layers of different materials
````

Here the counting begins from the lower layer, so in [] the first one belongs to the top layer and the second one to the lower layer.

Then you can choose which calculation method you want to use for layers. For example, sometimes you may use SLA method for a top layer to reduce the calculation time, so it means that it's possible to use different methods to calculate different layers.
```
Methods = {'SLA'};
Methods = {'SLA','NS'};
Methods = {'NS','NS'};
```
	
NS method is the default one. The number of pointed methods should coincide with the number of layers.

Now we can start the MultiLayer class:
```
Rml = ReflectionMultiLayer(Layers,Methods); %the reflection case
Rml = TransmissionMultiLayer(Layers,Methods); %the transmission case

Rml.energy_mesh_full = Mat.E0-100:step:Mat.E0+10; %to set the energy mesh you want to calculate (if you don't need to calculate the entire spectrum)
Rml.sigma_gauss = 0.8; % to convolute with the gauss function
Rml.Calculate; % etc, you can find all other things in example codes
```

### XPS
	
Layers and Methods can be set in the same way. The only one thing you need to add here isthe shell you would like to calculate. This consists of the next structure:
```
MS.Au={'4s1/2'};
MS.Si={'2s1/2'};
```
For example, we have a layer of Au on the top of a Si layer. So we need to say that we would like to calculate 4s1/2 Au and 2s1/2 Si.

Then,
```
Rml = PESMultiLayer(Layers,MS,Methods);
```
## To calculate additional DECS data with ELSEPA code

The source paper https://www.sciencedirect.com/science/article/pii/S0010465504004795?via%3Dihub
The class /src/MaterialDatabase/Data/ElsepaRunner.m is used to run the calculation of DECS data.

For each material, elastic (and not only) data need to be calculated for, a file Make_<nameofmaterial>.m must be created.
See examples in the folder /src/MaterialDatabase/Data/.
	
Let's take as an example Make_Au.m.
The line 
```
E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
```
determines the range of incident energies for which elastic scattering data are gonna be calculated.
And
```
[data] = ElsepaRunner.RunElsepa(Au.Z,E0);
```
runs the ELSEPA code for a certain atomic number Z and all incident energies E0.

To update the database MaterialData.mat MaterialData_script.m must be run with Make_<nameofmaterial>.m files included.
