import pappt
from unittest import TestSuite, findTestCases, TextTestRunner
test_suite = TestSuite()
def addTestCasesFromModule(module):
    test_suite.addTests(findTestCases(module))

addTestCasesFromModule(pappt.kvextract)
addTestCasesFromModule(pappt.lagrangian)
addTestCasesFromModule(pappt.pappt)

def run():
	TextTestRunner(verbosity=2).run(test_suite)

if __name__=='__main__':
	run()

