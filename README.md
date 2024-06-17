# Specific bond reaction 

Kimmdy reaction plugin to form specific bonds. Bonds can only occur between reactive particles (specified in the plumed.dat input file), which are not already forming a specific bond, and which are within a specified cutoff distance.

Fork from: reaction template for KIMMDY

## Installation
Should get installed together with kimmdy. If you want to install it separatly: 

* Download
* `pip install -e ./`  or for development: `pip install -r requirements.txt`

## Making your own
* Implement your reaction as a subclass of `kimmdy.reaction.Reaction`
* Register your Reaction class in the  **[options.entry_points]** section in the setup.cfg. The name you give here must match the entry in the config.yml for Kimmdy!




