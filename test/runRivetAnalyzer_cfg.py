import FWCore.ParameterSet.Config as cms

process = cms.Process("runRivetAnalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FE35E100-7544-E311-8869-7845C4FC36AD.root',
        #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FE608E15-7246-E311-8197-00266CF9AB88.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FEA8DB51-F343-E311-A964-848F69FD2484.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FEAD2571-1144-E311-B217-00266CF279F8.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FEB69B4F-4F44-E311-AE3B-00A0D1EEF4F8.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FEC65AA8-E643-E311-8E83-7845C4FC3CA1.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FECF8EA4-3E44-E311-848E-00A0D1EE8D00.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FECFC705-C046-E311-B9AF-008CFA001F78.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FED1F9F4-B744-E311-AFD6-00266CFAEA68.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FEEA3B27-EA43-E311-8DC2-00266CF24EEC.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/0037BFA7-D943-E311-8FA3-00266CF9C018.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/003CFA73-E945-E311-8917-00266CFAE318.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/005B4D11-F543-E311-B865-008CFA010D18.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/00817B05-8D45-E311-885C-848F69FD29DF.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/00903D2F-3E44-E311-8AAB-00266CF9B274.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/009499C5-CF43-E311-9C8E-7845C4FC3647.root', 
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/00A00D59-A244-E311-AAA3-00266CFAE810.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/00A610BE-4D44-E311-BAAD-848F69FD4667.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/00CB37BF-2745-E311-BF45-008CFA008D0C.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/00E2062C-9E44-E311-9EA5-00266CF253C4.root',
                                ### START m0-700 ###
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0039166F-C9CA-E211-9B10-0026189438E0.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/003AFE90-18CB-E211-AC56-002618943980.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/004D8BE1-AFCA-E211-A6F2-0026189438B9.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/005ABB25-FCCA-E211-BD1B-003048678E6E.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/007C2BC5-C5CA-E211-B4D1-003048679188.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/009A88E3-D6CA-E211-9220-003048678BB2.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/00B012B7-7CCA-E211-B974-003048678E52.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/00C1D3A3-BDCA-E211-8365-00261894394D.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/00EB5C96-D1CA-E211-99C1-0025905964BC.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/020CF657-DBCA-E211-AA85-00261894383B.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/021F51AF-FCCA-E211-9453-0026189438CB.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0229107C-A8CA-E211-AFCF-003048678BAC.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/023443D1-FBCA-E211-9E31-0026189438F5.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/02577C0E-C2CA-E211-8E6C-002618943906.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/027A4C33-B9CA-E211-88F3-0030486790A6.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/028C3515-0ECB-E211-9EDD-0030486791DC.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/028E2EF8-CBCA-E211-8343-003048FFCBB0.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/02B523DF-EBCA-E211-A597-0026189438F5.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/02D7EB86-AFCA-E211-9EB9-0026189438AF.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/02D95E21-E6CA-E211-9261-00261894386E.root',

                                ### START m700-1000 ###
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/002AF2F0-55E5-E211-ACA9-20CF3027A5B0.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/006FB40D-87E5-E211-B57E-001EC9D7F67B.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/02D1016B-7DE5-E211-BB84-20CF300E9EAD.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/042D97C0-C0E5-E211-95AC-00259073E4B6.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/086277F6-7BE5-E211-AF4B-485B39800B67.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/08A67227-BBE5-E211-B568-20CF3027A5FE.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0AAAB8E7-B3E5-E211-A8CC-E0CB4E29C4FF.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0C0C8B40-7AE5-E211-89E6-00261834B54C.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0C762741-73E5-E211-ABD9-001EC9D7F207.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0E3BD284-68E5-E211-BB2E-20CF3027A5DC.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0EFF2492-80E5-E211-BE99-20CF305B0582.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/1004E4CA-90E5-E211-A41F-485B39800B95.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/1028D409-B6E5-E211-8E08-E0CB4EA0A8FA.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/103591CC-BAE5-E211-8391-001EC9D8B986.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/104CA20A-9BE5-E211-B63C-20CF3027A5D1.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/1065398E-B7E5-E211-B7EE-20CF300E9ECB.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/128309BA-B4E5-E211-8B54-E0CB4EA0A8FA.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/129D33F3-B4E5-E211-A6B9-20CF3019DEE8.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/12B018D3-7DE5-E211-9B6C-00259073E474.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/12EC144B-5AE5-E211-9090-E0CB4E19F9A9.root',

                                ### m1000-Inf
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/000AF37A-7CE5-E211-A12F-00259073E31C.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/02A9CC2C-B4E5-E211-8DDE-90E6BA19A227.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0A2BB924-B4E5-E211-ADCD-20CF3027A564.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0C23B8C4-B6E5-E211-A605-00259073E3D0.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0C3D6A1F-6DE5-E211-9D5C-00261834B53C.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0C695C72-6FE5-E211-8435-90E6BAE8CC13.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/0C94CA53-D6E5-E211-BDB9-00259074AEDE.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/1045F5DD-D0E5-E211-ACDD-00259073E4DA.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/128ED289-B6E5-E211-B375-20CF3019DEE8.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/1606147D-88E5-E211-904A-001EC9D8D091.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/163F07E4-B3E5-E211-80F1-20CF305B057A.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/1835808E-75E5-E211-BC33-00259073E51A.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/1E1B78B4-BBE5-E211-BD53-20CF3027A5FE.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/1EE7CC6D-DCE5-E211-92AB-00259073E34A.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/20D0763D-B4E5-E211-A8F5-E0CB4EA0A8FA.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/26A20921-B6E5-E211-B6B0-0025907750A0.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/26EEDCAB-B6E5-E211-8170-20CF305616E0.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/28B3FEB0-7DE5-E211-A597-00259073E31C.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/28C3AC79-B5E5-E211-A427-90E6BA19A227.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_AUET2_8TeV-powheg-herwig/AODSIM/PU_S10_START53_V19-v1/10000/2AD1A7FF-B5E5-E211-BD7A-00259073E42E.root',

                                ### POWHEG + PYTHIA
                                ### m0-700
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/00052A2F-9901-E211-A26E-001A928116D6.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/001467F0-8601-E211-AB04-0026189438C2.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/0060FD65-C501-E211-A46F-001A92971B82.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/0095DA1D-B001-E211-A5F3-003048FFCC1E.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/009EC8EA-BE01-E211-A51F-002618943972.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/00B1850A-5C01-E211-B439-00261894383B.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/00E84574-AF01-E211-A72B-0026189437EC.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/023AFF7E-7901-E211-AB08-003048FFD720.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/024B7297-6401-E211-AF0E-002618943826.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/025D367D-7601-E211-84C9-001A9281173C.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/02684E78-7601-E211-B54E-001BFCDBD15E.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/027C04EA-A801-E211-926A-003048FFD71A.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/0290D012-9601-E211-BACF-001A928116E8.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/02996621-8901-E211-B3A1-003048679150.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/02A9FFF0-B501-E211-8620-002354EF3BD2.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/02CB5D8A-9701-E211-8273-001A92971B8C.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/02EB6936-6501-E211-8FED-002618943980.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/02FC5D64-5501-E211-BE5C-00261894380B.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/040CC435-8401-E211-8077-003048FFCB6A.root',
                                'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/04128E9C-B601-E211-9AB9-003048679266.root',
                                ### m700to1000
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v1/00000/005D3E2B-2909-E211-94DC-001E67397D7D.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-700to1000_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v1/00000/00C8600F-4E09-E211-ABEF-001E67397CAB.root',
                                ### m1000toInf
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v1/00000/00265667-E408-E211-A1C7-0030487F1F23.root',
                                #'root://cmsxrootd.fnal.gov//store/mc/Summer12_DR53X/TT_Mtt-1000toInf_CT10_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v1/00000/024B716A-E508-E211-892F-0030487D5E9D.root',

                            )
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("genParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
)
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2015_I1388555')
process.rivetAnalyzer.OutputFile = "MC_PowhegPythia_mInc.yoda"

process.p = cms.Path(process.generator*process.rivetAnalyzer)


