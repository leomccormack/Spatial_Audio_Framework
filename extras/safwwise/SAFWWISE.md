# safwwise

**extras/safwwise** is an initiative tasked with demonstrating how SAF may be integrated into [Wwise](https://www.audiokinetic.com/products/wwise/) plug-ins. The work should be viewed very much as a "proof-of-concept", and it is far from a finished product. However, it has been released as an open-source project, in the hope that it may prove useful for those interested in bringing SAF functionality to Wwise. 

## plugin_template

A tutorial which shows how to set up Wwise SDK and the required build tools to build a Wwise plugin for Windows. The plugin can be used in Wwise authoring or in Wwise-Unity integration (the latter isn't covered fully in depth but should be enough to get started if you know Unity and Wwise).

The actual plugin doesn't do anything but can be used as a template or a "hello world" type of test to learn the development tool environment and the workflow. This is meant to get a plugin built and up and running quickly. There doesn't seem to be a lot of material covering the workflow step-by-step so hopefully this will be useful for someone. For actual development and customization to your use case and platform, please refer to Wwise SDK documentation:

More information can be found [here](https://github.com/SAFrelated/plugin_template).

##  ambiBIN_generic

Another tutorial (which is quite experimental). It currently does not fully conform to Wwise requirements and is Windows 10 specific, but it may provide a reasonable starting point for tests and experiments with the Spatial Audio Framework (SAF) integrated within a Wwise effect plugin.

This project builds a somewhat working Authoring plugin which can be tested in Wwise, but is not yet fully functional in Wwise Unity Integration (shared plugin). The Authoring plugin has a generic Wwise GUI and the DSP functionality is similar to the SPARTA suite plugin AmbiBIN (without SOFA or OSC support).

More information can be found [here](https://github.com/SAFrelated/ambiBIN_generic).

## Contributing

Contributions are very much welcomed and encouraged!
