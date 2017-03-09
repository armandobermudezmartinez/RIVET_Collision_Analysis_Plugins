# Basic workflow

* Fork the repository
* Add a preliminary version of your plugin
* Create a work-in-progress (WIP) merge request
* Receive comments and finalize plugin
* Add validation plots to the merge request
* If code and validation is ok, CMS Rivet contacts will send the plugin to the Rivet authors

# Naming convention

* Papers: `CMS_{YYYY}_I{inspire-id}.cc`
* PAS: `CMS_{YYYY}_PAS_TOP_{YY}_{NNN}.cc`

If you do not have an inspire-id for your paper yet, please use the PAS naming convention.

If your PAS has an inspire-id, don't use it. Use the PAS naming convention.

# Plot style

Add `Title=CMS, xx TeV, process [<extra information>]` to your plots to make the CMS origin easily identifiable if theorists use the plugin.

# Coding style

Please follow the Rivet coding style documented here:
* https://rivet.hepforge.org/trac/wiki/WritingAnAnalysis
* https://rivet.hepforge.org/trac/wiki/CodingStyleGuide
