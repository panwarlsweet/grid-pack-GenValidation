import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/work/l/lata/MC/HIG-RunIIFall17wmLHEGS-00429.root'
                )
                            )



ntuple_genHiggs = cms.PSet(
     NtupleName = cms.string('NtupleGenJet'),
     genParticles = cms.InputTag('genParticles'),
)
ntuple_genjet = cms.PSet(
    NtupleName = cms.string('NtupleGenJet'),
    GenJets = cms.InputTag('ak4GenJets')
)
process.demo = cms.EDAnalyzer(
    "NtupleGenJet",
    Ntuples = cms.VPSet(
	ntuple_genHiggs,
	ntuple_genjet,
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_tree_bbH_ybyt_genpar.root"
))
process.p = cms.Path(process.demo)
process.MessageLogger.cerr.FwkReport.reportEvery = 100


