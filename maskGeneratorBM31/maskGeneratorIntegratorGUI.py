# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'guiLayout.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.
if __name__ == '__main__':
    import maskGeneratorIntegraterCdTe
    import maskGeneratorCdTe_recursive
else:
    from . import maskGeneratorIntegraterCdTe
    from . import maskGeneratorCdTe_recursive
from PyQt6 import QtCore, QtGui, QtWidgets
import os,time 

class Worker(QtCore.QThread):
    def __init__(self,direc,poni, mask,gainFile,recursePattern):
        super(Worker,self).__init__()
        self.direc = direc
        self.poni = poni
        self.mask = mask
        self.gainFile = gainFile
        self.recursePattern = recursePattern
    def run(self):
        self.running = True
        fileList = []
        while self.running:
            try:
                fileList, runningFull = maskGeneratorCdTe_recursive.run(self.direc,self.direc,self.poni,self.mask, self.gainFile, 
                                                                    self.recursePattern, fileList)
            except OSError as e:
                print(e)
                print('stopping')
                self.stop()
                return
            time.sleep(1)
            if runningFull:
                print('looking for new files')
    def stop(self):
        self.running = False

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(625, 331)
        MainWindow.setWindowTitle("imagin (Image MAsk Generator INtegrator)")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")

        self.configFile = f'{os.path.dirname(os.path.realpath(__file__))}/guiConfig.log'

        self.directoryBox = QtWidgets.QLineEdit(self.centralwidget)
        self.directoryBox.setGeometry(QtCore.QRect(10, 20, 431, 20))
        self.directoryBox.setObjectName("directoryBox")

        self.dirButton = QtWidgets.QPushButton(self.centralwidget)
        self.dirButton.setGeometry(QtCore.QRect(450, 20, 21, 21))
        self.dirButton.setObjectName("dirButton")
        self.dirButton.setText("...")

        self.poniBox = QtWidgets.QLineEdit(self.centralwidget)
        self.poniBox.setGeometry(QtCore.QRect(10, 60, 431, 20))
        self.poniBox.setObjectName("poniBox")

        self.maskBox = QtWidgets.QLineEdit(self.centralwidget)
        self.maskBox.setGeometry(QtCore.QRect(10, 100, 431, 20))
        self.maskBox.setObjectName("maskBox")

        self.gainMapBox = QtWidgets.QLineEdit(self.centralwidget)
        self.gainMapBox.setGeometry(QtCore.QRect(10, 140, 431, 20))
        self.gainMapBox.setObjectName("gainMapBox")

        self.runButton = QtWidgets.QPushButton(self.centralwidget)
        self.runButton.setGeometry(QtCore.QRect(340, 230, 75, 23))
        self.runButton.setObjectName("runButton")
        self.runButton.setText("Run")

        self.stopButton = QtWidgets.QPushButton(self.centralwidget)
        self.stopButton.setGeometry(QtCore.QRect(420, 230, 75, 23))
        self.stopButton.setObjectName("stopButton")
        self.stopButton.setText("Stop")
        self.stopButton.setEnabled(False)

        self.maskButton = QtWidgets.QPushButton(self.centralwidget)
        self.maskButton.setGeometry(QtCore.QRect(450, 100, 21, 21))
        self.maskButton.setObjectName("maskButton")
        self.maskButton.setText("...")

        self.poniButton = QtWidgets.QPushButton(self.centralwidget)
        self.poniButton.setGeometry(QtCore.QRect(450, 60, 21, 21))
        self.poniButton.setObjectName("poniButton")
        self.poniButton.setText("...")
        
        self.gainMapButton = QtWidgets.QPushButton(self.centralwidget)
        self.gainMapButton.setGeometry(QtCore.QRect(450, 140, 21, 21))
        self.gainMapButton.setObjectName("gainMapButton")
        self.gainMapButton.setText("...")

        self.dirLabel = QtWidgets.QLabel(self.centralwidget)
        self.dirLabel.setGeometry(QtCore.QRect(480, 20, 47, 16))
        self.dirLabel.setObjectName("dirLabel")
        self.dirLabel.setText( "directory")

        self.poniLabel = QtWidgets.QLabel(self.centralwidget)
        self.poniLabel.setGeometry(QtCore.QRect(480, 60, 47, 16))
        self.poniLabel.setObjectName("poniLabel")
        self.poniLabel.setText( "poni file")

        self.maskLabel = QtWidgets.QLabel(self.centralwidget)
        self.maskLabel.setGeometry(QtCore.QRect(480, 100, 47, 13))
        self.maskLabel.setObjectName("maskLabel")
        self.maskLabel.setText( "mask file")

        self.gainLabel = QtWidgets.QLabel(self.centralwidget)
        self.gainLabel.setGeometry(QtCore.QRect(480, 140, 71, 16))
        self.gainLabel.setObjectName("gainLabel")
        self.gainLabel.setText( "gain map file")   

        self.recurseBox = QtWidgets.QCheckBox(self.centralwidget)
        self.recurseBox.setGeometry(QtCore.QRect(20, 180, 20, 20))
        self.recurseBox.setObjectName("recurseBox")
        self.recurseBox.setText("run recursively\nin loop")
        self.recurseBox.adjustSize()

        self.recursePatternBox = QtWidgets.QLineEdit(self.centralwidget)
        self.recursePatternBox.setGeometry(QtCore.QRect(140, 180, 120, 20))
        self.recursePatternBox.setObjectName("recursePatternBox") 

        self.recurseLabel = QtWidgets.QLabel(self.centralwidget)
        self.recurseLabel.setGeometry(QtCore.QRect(280, 170, 71, 16))
        self.recurseLabel.setObjectName("recurseLabel")
        self.recurseLabel.setText( "recurse path pattern (only searches\ndirectories containing this pattern)") 
        self.recurseLabel.adjustSize()

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 625, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)     
        self.updateConfig()
        if os.path.exists(self.configFile):
            self.readConfigFile() 

        self.dirButton.clicked.connect(self.selectFolder)
        self.poniButton.clicked.connect(self.selectPoni)
        self.maskButton.clicked.connect(self.selectMask)
        self.gainMapButton.clicked.connect(self.selectGainFile)
        self.runButton.clicked.connect(self.run)
        self.stopButton.clicked.connect(self.stopWorker)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate

    def selectFolder(self):
        if self.directoryBox.text():
            currentdir = self.directoryBox.text()
        else:
            currentdir = '.'
        dialog = QtWidgets.QFileDialog.getExistingDirectory(caption = 'select folder', directory=currentdir)
        if dialog:
            self.directoryBox.setText(dialog)
            self.updateConfigFile()
    def selectPoni(self):
        if self.poniBox.text():
            currentdir = os.path.dirname(self.poniBox.text())
        else:
            currentdir = '.'
        filter = "data file (*.poni)"
        dialog = QtWidgets.QFileDialog.getOpenFileName(caption = 'select poni file',filter = filter, directory=currentdir)
        if dialog[0]:
            self.poniBox.setText(dialog[0])
            self.updateConfigFile()

    def selectMask(self):
        if self.maskBox.text():
            currentdir = os.path.dirname(self.maskBox.text())
        else:
            currentdir = '.'
        filter = "data file (*.edf)"
        dialog = QtWidgets.QFileDialog.getOpenFileName(caption = 'select mask file',	filter = filter, directory=currentdir)
        if dialog[0]:
            self.maskBox.setText(dialog[0])
            self.updateConfigFile()

    def selectGainFile(self):
        if self.gainMapBox.text():
            currentdir = os.path.dirname(self.gainMapBox.text())
        else:
            currentdir = '.'
        print(currentdir)
        filter = "data file (*.edf)"
        dialog = QtWidgets.QFileDialog.getOpenFileName(caption = 'select gain map file',	filter = filter, directory=currentdir)
        if dialog[0]:
            self.gainMapBox.setText(dialog[0])
            self.updateConfig()
            self.updateConfigFile()

    def updateConfig(self):
        self.configDct = {self.directoryBox.objectName() :[ self.directoryBox , self.directoryBox.text()],
                     self.poniBox.objectName(): [self.poniBox, self.poniBox.text()],
                     self.maskBox.objectName(): [self.maskBox, self.maskBox.text()],
                     self.gainMapBox.objectName(): [self.gainMapBox,self.gainMapBox.text()],
                     self.recursePatternBox.objectName(): [self.recursePatternBox, self.recursePatternBox.text()]}
        
    def updateConfigFile(self):
        self.updateConfig()
        string = ''
        for item in self.configDct:
            string += f'{item};{self.configDct[item][1]}\n'
        f = open(self.configFile,'w')
        f.write(string)
        f.close()

    def readConfigFile(self):
        f=open(self.configFile,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            line = line.replace('\n','')
            if not line:
                continue
            linesplit = line.split(';')
            self.configDct[linesplit[0]][0].setText(linesplit[1])
        self.updateConfig()
    
    def run(self):
        self.updateConfigFile()
        print('running')
        direc = self.directoryBox.text()
        dest = direc
        poni = self.poniBox.text()
        mask = self.maskBox.text()
        gainFile = self.gainMapBox.text()
        stringDct = {poni:'poni file',mask:'mask file', gainFile: 'gain file', direc: 'directory'}
        pars = [direc,poni,mask,gainFile]
        for par in pars:
            if not par:
                print(f'input {stringDct[par]}')
                return
            elif not os.path.exists(par):
                print(f'{stringDct[par]} not found')
                return

        try:
            if self.recurseBox.isChecked():
                dirpattern = self.recursePatternBox.text()
                self.startWorker(direc,poni,mask,gainFile,dirpattern)
                #maskGeneratorCdTe_recursive.run(direc,dest,poni,mask,gainFile, dirpattern)
            else:
                maskGeneratorIntegraterCdTe.run(direc,dest,poni,mask,gainFile)
            print('finished')
        except IndexError:
            print('no valid files')
            return
        except FileNotFoundError as e:
            print(e)
            return
        
    def startWorker(self, direc, poni,mask,gainFile,recursePattern):
        self.runButton.setEnabled(False)
        self.stopButton.setEnabled(True)
        self.thread = Worker(direc, poni,mask,gainFile,recursePattern)
        self.thread.start()

    def stopWorker(self):
        self.thread.stop()
        self.thread.quit()
        self.thread.wait()
        self.runButton.setEnabled(True)
        self.stopButton.setEnabled(False)
        print('stopping')

            
def main():
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
