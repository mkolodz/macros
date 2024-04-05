#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TArrow.h"

#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TString.h"
#include "TMatrixD.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h> 

#define NBINS 20
#define NBINSLESS 10
void readFromSQLAndPlotResultsForFiberStacks()
{
        gStyle->SetPalette(kBird);
//     gStyle->SetOptStat(0);
        
    sqlite3 *db;
    sqlite3_stmt *stmt;
    char *zErrMsg = 0;
    const char *sql = "select ATT_COMB,ATT_COMB_ERR from ATTENUATION_LENGTH";

    
    TH1D * hAttComb = new TH1D("hAttComb", "hAttComb", NBINS, 50, 600);    
    TH1D * hAttCombErr = new TH1D("hAttCombErr", "hAttCombErr", NBINS, 0, 25);  
    TH1D * hS0 = new TH1D("hS0", "hS0", NBINS, 100, 600);    
    TH1D * hS0Err = new TH1D("hS0Err", "hS0Err", NBINS, 15, 50);  
    TH1D * hLambda = new TH1D("hLambda", "hLambda", NBINS, 50, 200);    
    TH1D * hLambdaErr = new TH1D("hLambdaErr", "hLambdaErr", NBINS, 0, 20);  
    TH1D * hEtaR = new TH1D("hEtaR", "hEtaR", NBINS, 0, 25);    
    TH1D * hEtaRErr = new TH1D("hEtaRErr", "hEtaRErr", NBINS, 0, 15);  
    TH1D * hEtaL = new TH1D("hEtaL", "hEtaL", NBINS, 0, 35);    
    TH1D * hEtaLErr = new TH1D("hEtaLErr", "hEtaLErr", NBINS, 0, 15);  
    TH1D * hKsi = new TH1D("hKsi", "hKsi", NBINS, 0, 5);    
    TH1D * hKsiErr = new TH1D("hKsiErr", "hKsiErr", NBINS, 0, 0.3);  
    TH1D * hChi2NDF = new TH1D("hChi2NDF", "hChi2NDF", NBINS, 0, 3);  
    
    TH1D * hPositionResPol1 = new TH1D("hPositionResPol1", "hPositionResPol1", NBINS, 0, 120);    
    TH1D * hPositionResPol1Err = new TH1D("hPositionResPol1Err", "hPositionResPol1Err", NBINS, 0, 1.2);    
    TH1D * hPositionResPol1All = new TH1D("hPositionResPol1All", "hPositionResPol1All", NBINS, 0, 200);    
    TH1D * hPositionResPol1AllErr = new TH1D("hPositionResPol1AllErr", "hPositionResPol1AllErr", NBINS, 0, 1);    

    TH1D * hEnresAv = new TH1D("hEnresAv", "hEnresAv", NBINS, 5,20);    
    TH1D * hEnresAvErr = new TH1D("hEnresAvErr", "hEnresAvErr", NBINS, 0, 0.5);    
    TH1D * hEnresCh0 = new TH1D("hEnresCh0", "hEnresCh0", NBINS, 5,20);    
    TH1D * hEnresCh0Err = new TH1D("hEnresCh0Err", "hEnresCh0Err", NBINS,0, 0.5);    
    TH1D * hEnresCh1 = new TH1D("hEnresCh1", "hEnresCh1", NBINS, 5,20);    
    TH1D * hEnresCh1Err = new TH1D("hEnresCh1Err", "hEnresCh1Err", NBINS,0, 0.5); 
    
    TH1D * hAttCombESR = new TH1D("hAttCombESR", "hAttCombESR", NBINS, 50, 600);    
    TH1D * hAttCombErrESR = new TH1D("hAttCombErrESR", "hAttCombErrESR", NBINS, 0, 25);  
    TH1D * hS0ESR = new TH1D("hS0ESR", "hS0ESR", NBINS, 100, 600);    
    TH1D * hS0ErrESR = new TH1D("hS0ErrESR", "hS0ErrESR", NBINS, 15, 50);  
    TH1D * hLambdaESR = new TH1D("hLambdaESR", "hLambdaESR", NBINS, 50, 200);    
    TH1D * hLambdaErrESR = new TH1D("hLambdaErrESR", "hLambdaErrESR", NBINS, 0, 20);  
    TH1D * hEtaRESR = new TH1D("hEtaRESR", "hEtaRESR", NBINS, 0, 25);    
    TH1D * hEtaRErrESR = new TH1D("hEtaRErrESR", "hEtaRErrESR", NBINS, 0, 15);  
    TH1D * hEtaLESR = new TH1D("hEtaLESR", "hEtaLESR", NBINS, 0, 35);    
    TH1D * hEtaLErrESR = new TH1D("hEtaLErrESR", "hEtaLErrESR", NBINS, 0, 15);  
    TH1D * hKsiESR = new TH1D("hKsiESR", "hKsiESR", NBINS, 0, 5);    
    TH1D * hKsiErrESR = new TH1D("hKsiErrESR", "hKsiErrESR", NBINS, 0, 0.3); 
    TH1D * hChi2NDFESR = new TH1D("hChi2NDFESR", "hChi2NDFESR", NBINS, 0, 3);  
    
    TH1D * hPositionResPol1ESR = new TH1D("hPositionResPol1ESR", "hPositionResPol1ESR", NBINS, 0, 120);    
    TH1D * hPositionResPol1ErrESR = new TH1D("hPositionResPol1ErrESR", "hPositionResPol1ErrESR", NBINS, 0, 1.2);    
    TH1D * hPositionResPol1AllESR = new TH1D("hPositionResPol1AllESR", "hPositionResPol1AllESR", NBINS, 0, 200);    
    TH1D * hPositionResPol1AllErrESR = new TH1D("hPositionResPol1AllErrESR", "hPositionResPol1AllErrESR", NBINS, 0, 1);  

    TH1D * hEnresAvESR = new TH1D("hEnresAvESR", "hEnresAvESR", NBINS, 5,20);    
    TH1D * hEnresAvErrESR = new TH1D("hEnresAvErrESR", "hEnresAvErrESR", NBINS, 0, 0.5);    
    TH1D * hEnresCh0ESR = new TH1D("hEnresCh0ESR", "hEnresCh0ESR", NBINS, 5,20);    
    TH1D * hEnresCh0ErrESR = new TH1D("hEnresCh0ErrESR", "hEnresCh0ErrESR", NBINS,0, 0.5);    
    TH1D * hEnresCh1ESR = new TH1D("hEnresCh1ESR", "hEnresCh1ESR", NBINS, 5,20);    
    TH1D * hEnresCh1ErrESR = new TH1D("hEnresCh1ErrESR", "hEnresCh1ErrESR", NBINS,0, 0.5); 
    
    sqlite3_open("/scratch1/gccb/data/results/FiberStackAl/FiberStackAlResults.db", &db);

    
    int rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        printf("error: %s", sqlite3_errmsg(db));
        return;
    }
    sqlite3_bind_double(stmt, 1, 1.); //bind values to prepared statements
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        double att_len = sqlite3_column_double(stmt, 0);
        hAttComb->Fill(att_len);
        double att_len_err = sqlite3_column_double(stmt, 1);
        hAttCombErr->Fill(att_len_err);
    }
    if (rc != SQLITE_DONE) {
        printf("error: %s", sqlite3_errmsg(db));
    }
    sqlite3_finalize(stmt);
    
    sql = "select S0,S0_ERR, LAMBDA, LAMBDA_ERR, ETAR, ETAR_ERR, ETAL, ETAL_ERR, KSI, KSI_ERR, CHI2NDF from  ATTENUATION_MODEL ";
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        printf("error: %s", sqlite3_errmsg(db));
        return;
    }
    sqlite3_bind_double(stmt, 1, 1.); //bind values to prepared statements
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        hS0->Fill(sqlite3_column_double(stmt, 0));
        hS0Err->Fill(sqlite3_column_double(stmt, 1));
        hLambda->Fill(sqlite3_column_double(stmt, 2));
        hLambdaErr->Fill(sqlite3_column_double(stmt, 3));
        hEtaR->Fill(sqlite3_column_double(stmt, 4));
        hEtaRErr->Fill(sqlite3_column_double(stmt, 5));
        hEtaL->Fill(sqlite3_column_double(stmt, 6));
        hEtaLErr->Fill(sqlite3_column_double(stmt, 7));
        hKsi->Fill(sqlite3_column_double(stmt, 8));
        hKsiErr->Fill(sqlite3_column_double(stmt, 9));
        hChi2NDF->Fill(sqlite3_column_double(stmt, 10));
    }
    if (rc != SQLITE_DONE) {
        printf("error: %s", sqlite3_errmsg(db));
    }
    sqlite3_finalize(stmt);
    
    sql="select POSITION_RES_POL1, POSITION_RES_POL1_ERR, POSITION_RES_POL1_ALL, POSITION_RES_POL1_ALL_ERR from  POSITION_RESOLUTION";
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        printf("error: %s", sqlite3_errmsg(db));
        return;
    }
    sqlite3_bind_double(stmt, 1, 1.); //bind values to prepared statements
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        hPositionResPol1->Fill(sqlite3_column_double(stmt, 0));
        hPositionResPol1Err->Fill(sqlite3_column_double(stmt, 1));
        hPositionResPol1All->Fill(sqlite3_column_double(stmt, 2));
        hPositionResPol1AllErr->Fill(sqlite3_column_double(stmt, 3));
    }
    if (rc != SQLITE_DONE) {
        printf("error: %s", sqlite3_errmsg(db));
    }
    sqlite3_finalize(stmt);
    
    sql="select ENRES_AV, ENRES_AV_ERR, ENRES_CH0, ENRES_CH0_ERR, ENRES_CH1, ENRES_CH1_ERR from  ENERGY_RESOLUTION";
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        printf("error: %s", sqlite3_errmsg(db));
        return;
    }
    sqlite3_bind_double(stmt, 1, 1.); //bind values to prepared statements
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        hEnresAv->Fill(sqlite3_column_double(stmt, 0));
        hEnresAvErr->Fill(sqlite3_column_double(stmt, 1));
        hEnresCh0->Fill(sqlite3_column_double(stmt, 2));
        hEnresCh0Err->Fill(sqlite3_column_double(stmt, 3));
        hEnresCh1->Fill(sqlite3_column_double(stmt, 4));
        hEnresCh1Err->Fill(sqlite3_column_double(stmt, 5));
    }
    if (rc != SQLITE_DONE) {
        printf("error: %s", sqlite3_errmsg(db));
    }
    sqlite3_finalize(stmt);
    
    sqlite3_close(db);
    
    
    sqlite3_open("/scratch1/gccb/data/results/FiberStackESR/FiberStackESRResults.db", &db);
    sql = "select ATT_COMB,ATT_COMB_ERR from ATTENUATION_LENGTH";
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        printf("error: %s", sqlite3_errmsg(db));
        return;
    }
    sqlite3_bind_double(stmt, 1, 1.); //bind values to prepared statements
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        double att_len = sqlite3_column_double(stmt, 0);
        hAttCombESR->Fill(att_len);
        double att_len_err = sqlite3_column_double(stmt, 1);
        hAttCombErrESR->Fill(att_len_err);
    }
    if (rc != SQLITE_DONE) {
        printf("error: %s", sqlite3_errmsg(db));
    }
    sqlite3_finalize(stmt);
    
    sql = "select S0,S0_ERR, LAMBDA, LAMBDA_ERR, ETAR, ETAR_ERR, ETAL, ETAL_ERR, KSI, KSI_ERR, CHI2NDF from  ATTENUATION_MODEL";
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        printf("error: %s", sqlite3_errmsg(db));
        return;
    }
    sqlite3_bind_double(stmt, 1, 1.); //bind values to prepared statements
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        hS0ESR->Fill(sqlite3_column_double(stmt, 0));
        hS0ErrESR->Fill(sqlite3_column_double(stmt, 1));
        hLambdaESR->Fill(sqlite3_column_double(stmt, 2));
        hLambdaErrESR->Fill(sqlite3_column_double(stmt, 3));
        hEtaRESR->Fill(sqlite3_column_double(stmt, 4));
        hEtaRErrESR->Fill(sqlite3_column_double(stmt, 5));
        hEtaLESR->Fill(sqlite3_column_double(stmt, 6));
        hEtaLErrESR->Fill(sqlite3_column_double(stmt, 7));
        hKsiESR->Fill(sqlite3_column_double(stmt, 8));
        hKsiErrESR->Fill(sqlite3_column_double(stmt, 9));
        hChi2NDFESR->Fill(sqlite3_column_double(stmt, 10));
    }
    if (rc != SQLITE_DONE) {
        printf("error: %s", sqlite3_errmsg(db));
    }
    sqlite3_finalize(stmt);
    
        sql="select POSITION_RES_POL1, POSITION_RES_POL1_ERR, POSITION_RES_POL1_ALL, POSITION_RES_POL1_ALL_ERR from  POSITION_RESOLUTION";
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        printf("error: %s", sqlite3_errmsg(db));
        return;
    }
    sqlite3_bind_double(stmt, 1, 1.); //bind values to prepared statements
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        hPositionResPol1ESR->Fill(sqlite3_column_double(stmt, 0));
        hPositionResPol1ErrESR->Fill(sqlite3_column_double(stmt, 1));
        hPositionResPol1AllESR->Fill(sqlite3_column_double(stmt, 2));
        hPositionResPol1AllErrESR->Fill(sqlite3_column_double(stmt, 3));
    }
    if (rc != SQLITE_DONE) {
        printf("error: %s", sqlite3_errmsg(db));
    }
    sqlite3_finalize(stmt);
    
    sql="select ENRES_AV, ENRES_AV_ERR, ENRES_CH0, ENRES_CH0_ERR, ENRES_CH1, ENRES_CH1_ERR from  ENERGY_RESOLUTION";
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        printf("error: %s", sqlite3_errmsg(db));
        return;
    }
    sqlite3_bind_double(stmt, 1, 1.); //bind values to prepared statements
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
        hEnresAvESR->Fill(sqlite3_column_double(stmt, 0));
        hEnresAvErrESR->Fill(sqlite3_column_double(stmt, 1));
        hEnresCh0ESR->Fill(sqlite3_column_double(stmt, 2));
        hEnresCh0ErrESR->Fill(sqlite3_column_double(stmt, 3));
        hEnresCh1ESR->Fill(sqlite3_column_double(stmt, 4));
        hEnresCh1ErrESR->Fill(sqlite3_column_double(stmt, 5));
    }
    if (rc != SQLITE_DONE) {
        printf("error: %s", sqlite3_errmsg(db));
    }
    sqlite3_finalize(stmt);
    
    
    sqlite3_close(db);
    
    auto canAtt = new TCanvas("Att","Att",2000,1500);
    
    canAtt->SetLeftMargin(0.4);
    canAtt->Divide(4,4);
    
    canAtt->cd(1);
    gPad->SetLeftMargin(0.15);
    hAttComb->Draw();
    hAttCombESR->Draw("same");
    hAttCombESR->SetLineColor(kRed);
    hAttComb->SetTitle("AttComb");
    hAttComb->GetYaxis()->SetRangeUser(0, hAttComb->GetMaximum() > hAttCombESR->GetMaximum() ? hAttComb->GetMaximum() + 1 : hAttCombESR->GetMaximum() + 1);

    
    canAtt->cd(2);
    gPad->SetLeftMargin(0.15);
    hAttCombErr->Draw();
    hAttCombErrESR->Draw("same");
    hAttCombErrESR->SetLineColor(kRed);
    hAttCombErr->SetTitle("AttCombErr");
    
    canAtt->cd(3);
    gPad->SetLeftMargin(0.15);
    hS0->Draw();
    hS0ESR->Draw("same");
    hS0ESR->SetLineColor(kRed);
    hS0->SetTitle("S0");
    
    canAtt->cd(4);
    gPad->SetLeftMargin(0.15);
    hS0Err->Draw();
    hS0ErrESR->Draw("same");
    hS0ErrESR->SetLineColor(kRed);
    hS0Err->SetTitle("S0Err");
    hS0Err->GetYaxis()->SetRangeUser(0, hS0Err->GetMaximum() > hS0ErrESR->GetMaximum() ? hS0Err->GetMaximum() + 1 : hS0ErrESR->GetMaximum() + 1);

        
    canAtt->cd(5);
    gPad->SetLeftMargin(0.15);
    hLambda->Draw();
    hLambdaESR->Draw("same");
    hLambdaESR->SetLineColor(kRed);
    hLambda->SetTitle("Lambda");
    hLambda->GetYaxis()->SetRangeUser(0, hLambda->GetMaximum() > hLambdaESR->GetMaximum() ? hLambda->GetMaximum() + 1 : hLambdaESR->GetMaximum() + 1);

    canAtt->cd(6);
    gPad->SetLeftMargin(0.15);
    hLambdaErr->Draw();
    hLambdaErrESR->Draw("same");
    hLambdaErrESR->SetLineColor(kRed);
    hLambdaErr->SetTitle("LambdaErr");
    hLambdaErr->GetYaxis()->SetRangeUser(0, hLambdaErr->GetMaximum() > hLambdaErrESR->GetMaximum() ? hLambdaErr->GetMaximum() + 1 : hLambdaErrESR->GetMaximum() + 1);

    canAtt->cd(7);
    gPad->SetLeftMargin(0.15);
    hEtaR->Draw();
    hEtaRESR->Draw("same");
    hEtaRESR->SetLineColor(kRed);
    hEtaR->SetTitle("EtaR");
    hEtaR->GetYaxis()->SetRangeUser(0, hEtaR->GetMaximum() > hEtaRESR->GetMaximum() ? hEtaR->GetMaximum() + 1 : hEtaRESR->GetMaximum() + 1);

    
    canAtt->cd(8);
    gPad->SetLeftMargin(0.15);
    hEtaRErr->Draw();
    hEtaRErrESR->Draw("same");
    hEtaRErrESR->SetLineColor(kRed);
    hEtaRErr->SetTitle("EtaRErr");
    hEtaRErr->GetYaxis()->SetRangeUser(0, hEtaRErr->GetMaximum() > hEtaRErrESR->GetMaximum() ? hEtaRErr->GetMaximum() + 1 : hEtaRErrESR->GetMaximum() + 1);

    
    canAtt->cd(9);
    gPad->SetLeftMargin(0.15);
    hEtaL->Draw();
    hEtaLESR->Draw("same");
    hEtaLESR->SetLineColor(kRed);
    hEtaL->SetTitle("EtaL");
    hEtaL->GetYaxis()->SetRangeUser(0, hEtaL->GetMaximum() > hEtaLESR->GetMaximum() ? hEtaL->GetMaximum() + 1 : hEtaLESR->GetMaximum() + 1);

    
    canAtt->cd(10);
    gPad->SetLeftMargin(0.15);
    hEtaLErr->Draw();
    hEtaLErrESR->Draw("same");
    hEtaLErrESR->SetLineColor(kRed);
    hEtaLErr->SetTitle("EtaLErr");
    hEtaLErr->GetYaxis()->SetRangeUser(0, hEtaLErr->GetMaximum() > hEtaLErrESR->GetMaximum() ? hEtaLErr->GetMaximum() + 1 : hEtaLErrESR->GetMaximum() + 1);

    
    canAtt->cd(11);
    gPad->SetLeftMargin(0.15);
    hKsi->Draw();
    hKsiESR->Draw("same");
    hKsiESR->SetLineColor(kRed);
    hKsi->SetTitle("Ksi");
    
    canAtt->cd(12);
    gPad->SetLeftMargin(0.15);
    hKsiErr->Draw();
    hKsiErrESR->Draw("same");
    hKsiErrESR->SetLineColor(kRed);
    hKsiErr->SetTitle("KsiErr");
    hKsiErr->GetYaxis()->SetRangeUser(0, hKsiErr->GetMaximum() > hKsiErrESR->GetMaximum() ? hKsiErr->GetMaximum() + 1 : hKsiErrESR->GetMaximum() + 1);

    canAtt->cd(13);
    gPad->SetLeftMargin(0.15);
    hChi2NDF->Draw();
    hChi2NDFESR->Draw("same");
    hChi2NDFESR->SetLineColor(kRed);
    hChi2NDF->SetTitle("Chi2NDF");
    hChi2NDF->GetYaxis()->SetRangeUser(0, hChi2NDF->GetMaximum() > hChi2NDFESR->GetMaximum() ? hChi2NDF->GetMaximum() + 1 : hChi2NDFESR->GetMaximum() + 1);


    auto canRes = new TCanvas("Res","Res",2000,1500);
    
    canRes->SetLeftMargin(0.4);
    canRes->Divide(4,3);

    canRes->cd(1);
    gPad->SetLeftMargin(0.15);
    hPositionResPol1->Draw();
    hPositionResPol1ESR->Draw("same");
    hPositionResPol1ESR->SetLineColor(kRed);
    hPositionResPol1->SetTitle("PositionResPol1");
    hPositionResPol1->GetYaxis()->SetRangeUser(0, hPositionResPol1->GetMaximum() > hPositionResPol1ESR->GetMaximum() ? hPositionResPol1->GetMaximum() + 1 : hPositionResPol1ESR->GetMaximum() + 1);
    
    canRes->cd(2);
    gPad->SetLeftMargin(0.15);
    hPositionResPol1Err->Draw();
    hPositionResPol1ErrESR->Draw("same");
    hPositionResPol1ErrESR->SetLineColor(kRed);
    hPositionResPol1Err->SetTitle("PositionResPolErr1");
    
        canRes->cd(3);
    gPad->SetLeftMargin(0.15);
    hPositionResPol1All->Draw();
    hPositionResPol1AllESR->Draw("same");
    hPositionResPol1AllESR->SetLineColor(kRed);
    hPositionResPol1All->SetTitle("PositionResPol1All");
    
        canRes->cd(4);
    gPad->SetLeftMargin(0.15);
    hPositionResPol1AllErr->Draw();
    hPositionResPol1AllErrESR->Draw("same");
    hPositionResPol1AllErrESR->SetLineColor(kRed);
    hPositionResPol1AllErr->SetTitle("PositionResPol1AllErr");
    
        canRes->cd(5);
    gPad->SetLeftMargin(0.15);
    hEnresAv->Draw();
    hEnresAvESR->Draw("same");
    hEnresAvESR->SetLineColor(kRed);
    hEnresAv->SetTitle("EnresAv");
        hEnresAv->GetYaxis()->SetRangeUser(0, hEnresAv->GetMaximum() > hEnresAvESR->GetMaximum() ? hEnresAv->GetMaximum() + 1 : hEnresAvESR->GetMaximum() + 1);

    
        canRes->cd(6);
    gPad->SetLeftMargin(0.15);
    hEnresAvErr->Draw();
    hEnresAvErrESR->Draw("same");
    hEnresAvErrESR->SetLineColor(kRed);
    hEnresAvErr->SetTitle("EnresAvErr");
        hEnresAvErr->GetYaxis()->SetRangeUser(0, hEnresAvErr->GetMaximum() > hEnresAvErrESR->GetMaximum() ? hEnresAvErr->GetMaximum() + 1 : hEnresAvErrESR->GetMaximum() + 1);

    
        canRes->cd(7);
    gPad->SetLeftMargin(0.15);
    hEnresCh0->Draw();
    hEnresCh0ESR->Draw("same");
    hEnresCh0ESR->SetLineColor(kRed);
    hEnresCh0->SetTitle("EnresCh0");
    
        canRes->cd(8);
    gPad->SetLeftMargin(0.15);
    hEnresCh0Err->Draw();
    hEnresCh0ErrESR->Draw("same");
    hEnresCh0ErrESR->SetLineColor(kRed);
    hEnresCh0Err->SetTitle("EnresCh0Err");
    
        canRes->cd(9);
    gPad->SetLeftMargin(0.15);
    hEnresCh1->Draw();
    hEnresCh1ESR->Draw("same");
    hEnresCh1ESR->SetLineColor(kRed);
    hEnresCh1->SetTitle("EnresCh1");
    
        canRes->cd(10);
    gPad->SetLeftMargin(0.15);
    hEnresCh1Err->Draw();
    hEnresCh1ErrESR->Draw("same");
    hEnresCh1ErrESR->SetLineColor(kRed);
    hEnresCh1Err->SetTitle("EnresCh1Err");

}
