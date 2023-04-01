/**
 * @file mainwindow.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-04-01
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "mainwindow.h"

#include <stdio.h>

#include <QtGui/QWheelEvent>
#include <QtCore/QRandomGenerator>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>

#include "ui_mainwindow.h"
#include "nodeandedge.h"

MainWindow::MainWindow(QWidget *parent):
    QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    theScene = new QGraphicsScene(ui->graphWidget);
    theScene->setItemIndexMethod(QGraphicsScene::NoIndex);
    theScene->setSceneRect(ui->graphWidget->geometry());
    ui->graphWidget->setScene(theScene);
    ui->graphWidget->setCacheMode(ui->graphWidget->CacheBackground);
    ui->graphWidget->setViewportUpdateMode(ui->graphWidget->BoundingRectViewportUpdate);
    ui->graphWidget->setRenderHint(QPainter::Antialiasing);
    ui->graphWidget->setTransformationAnchor(ui->graphWidget->AnchorUnderMouse);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_action_Open_triggered() {

}

void MainWindow::on_action_New_triggered() {

}

void MainWindow::on_action_Save_triggered() {
    if (filename == ""){
        auto fileUrl = QFileDialog::getSaveFileUrl();
        auto file = std::fopen(fileUrl.toLocalFile().toStdString().c_str(), "w");
        if (file == NULL) {
            QMessageBox::information(this, "Error", "Failed to open file!");
            return;
        }
        filename = fileUrl.toLocalFile().toStdString();
    }
    auto file = std::fopen(filename.c_str(), "w");
    save(file);
}

void MainWindow::on_actionSave_As_triggered() {
    auto fileUrl = QFileDialog::getSaveFileUrl();
    auto file = std::fopen(fileUrl.toLocalFile().toStdString().c_str(), "w");
    if (file == NULL) {
        QMessageBox::information(this, "Error", "Failed to open file!");
        return;
    }
    save(file);
}

void MainWindow::save(FILE* file) {
    const char __data_file_header[] =
        "# Modeled by designer. AUTO generated, DO NOT EDIT\n\n"            // header
        "%li\tatoms\n%li\tbonds\n\n1\tatom types\n1\tbond types\n\n"       // atom & bond number
        "0\t20\txlo xhi\n0\t20\tylo yhi\n0\t20\tzlo zhi\n\n"  // boundary
        "Masses\n\n1\t1\n";

    // header, check whether bond_relations2_.size() == 0 is for model3
    std::fprintf(file, __data_file_header, static_cast<size_t>(nodes.size()),
                 static_cast<size_t>(edges.size()));

    // Atoms, id type x y z
    std::fprintf(file, "\nAtoms\n\n");
    for (auto i : nodes)
        std::fprintf(file, "%s", i->print().c_str());
  
    // Bonds, id type a b
    std::fprintf(file, "\nBonds\n\n");
    for (auto i : edges)
        std::fprintf(file, "%s", i->print().c_str());
    std::fclose(file);
}

void MainWindow::on_actionAdd_Node_triggered() {
    Node *node = new Node(ui->graphWidget, nodes.size());
    nodes.append(node);
    node->setPos(50, 50);
    theScene->addItem(node);
}

void MainWindow::on_actionRemove_Node_triggered(){
    for (int i=nodes.size()-1; i>=0; i--)
    if (nodes[i]->selected) {
        for (int j=nodes[i]->edgeList.size()-1; j>=0; j--) {
            auto edge = nodes[i]->edgeList[j];
            theScene->removeItem(edge);
            edges.remove(edge->id);
            delete edge;
            for (int k=0; k<edges.size(); k++)
                edges[k]->id = k;
        }
        theScene->removeItem(nodes[i]);
        delete nodes[i];
        nodes.remove(i);
        for (int i=0; i<nodes.size(); i++) {
            nodes[i]->id = i;
        }
    }
}

Edge *MainWindow::findSelectedEdge(bool create) {
    Node *a = nullptr, *b = nullptr;
    Edge *edge = nullptr;

    // check whether just two nodes are selected
    for (int i=0; i<nodes.size(); i++) {
        if (nodes[i]->selected) {
            if (a == nullptr) {
                a = nodes[i];
                continue;
            }
            if (b == nullptr) {
                b = nodes[i];
                continue;
            }
            QMessageBox::information(this, "Error", "More than two nodes are selected!");
            return nullptr;
        }
    }
    if (a == nullptr || b == nullptr) {
        QMessageBox::information(this, "Error", "Less than two nodes are selected!");
        return nullptr;
    }

    // find whether node been created
    for (int i=0 ; i<edges.size(); i++) {
        if (edges[i]->isConnected(a) && edges[i]->isConnected(b))
            edge = edges[i];
    }

    if (create) {
        if (edge != nullptr) {
            QMessageBox::information(this, "Error", "Edge has been created!");
            return nullptr;
        }
        edge = new Edge(a, b, edges.size());
        return edge;
    }

    if (edge == nullptr)
        QMessageBox::information(this, "Error", "Edge hasn't been created!");
    return edge;
}

void MainWindow::on_actionAdd_Edge_triggered() {
    Edge* edge = findSelectedEdge(true);
    if (edge == nullptr)
        return;
    theScene->addItem(edge);
    edges.append(edge);
}

void MainWindow::on_actionRemove_Edge_triggered() {
    Edge* edge = findSelectedEdge(false);
    if (edge == nullptr)
        return;
    theScene->removeItem(edge);
    edges.remove(edge->id);
    delete edge;
    for (int i=0; i<edges.size(); i++)
        edges[i]->id = i;
}

