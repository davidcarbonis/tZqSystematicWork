import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100000)

process.load( "SimGeneral.HepPDTESSource.pythiapdt_cfi" )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ## Scale up files
        'file:/nfs/data/eepgadm/tZq/MC/2012/TTZJets_8TeV_scaleup/GEN/genFile0.root',
        'file:/nfs/data/eepgadm/tZq/MC/2012/TTZJets_8TeV_scaleup/GEN/genFile1.root',
        'file:/nfs/data/eepgadm/tZq/MC/2012/TTZJets_8TeV_scaleup/GEN/genFile2.root',

        ## Scale down files
        #'file:/nfs/data/eepgadm/tZq/MC/2012/TTZJets_8TeV_scaledown/GEN/genFile0.root',
        #'file:/nfs/data/eepgadm/tZq/MC/2012/TTZJets_8TeV_scaledown/GEN/genFile1.root',
        #'file:/nfs/data/eepgadm/tZq/MC/2012/TTZJets_8TeV_scaledown/GEN/genFile2.root',
    )
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("ScaleUpOutput.root")
)


process.analysis = cms.EDAnalyzer('TTZgenLevelAnalyser'
)


process.p = cms.Path(process.analysis)
