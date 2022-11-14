How to setup vercel.com:

* In vercel.com go create a new project.
* Import it from the git repository.
* Select personal account
* Give it a Project Name and apply these settings:
* * Framekwork: Other
* * Build Command: ./1_build.sh
* * Output Directory: ./export
* * Install Command: ./0_setup.sh
* * Root Directory: test/documentation/.vercel
* * * Enabled: Include source files outside of the Root Directory in the Build Step.
