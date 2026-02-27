import parmed as pmd

amber = pmd.load_file('disarcosine.top', 'disarcosine.crd')

amber.save('disarcosine_gmx.top')
amber.save('disarcosine_gmx.gro')
