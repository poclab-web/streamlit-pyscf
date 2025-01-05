

### style.css
streamlitのボタンの色の修正などをcustom cssを用いて行なった。

### solvents_epsilon.csv
solventsのPCMに設定するパラメーターのデータは、以下を参考にして行なった。
https://pyscf.org/user/solvent.html

Solvent parameters
The default solvent in the PCM module is water. When studying other types of solvents, you can consider to modify the dielectric parameter eps using the constants listed below:

import pyscf
mol = pyscf.M(atom='''
     C  0.    0.      -0.542
     O  0.    0.       0.677
     H  0.    0.935   -1.082
     H  0.   -0.935   -1.082''',
              basis='6-31g*', verbose=4)
mf = mol.RHF().PCM()
mf.with_solvent.eps = 32.613   # methanol
mf.run()

以下のページのSOLVENTSというタブから数値を得た。
https://gaussian.com/scrf/