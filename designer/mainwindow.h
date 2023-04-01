/**
 * @file mainwindow.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-04-01
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QGraphicsView>

class Node;
class Edge;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_action_Open_triggered();
    void on_action_New_triggered();
    void on_action_Save_triggered();
    void on_actionSave_As_triggered();
    void on_actionAdd_Node_triggered();
    void on_actionRemove_Node_triggered();
    void on_actionAdd_Edge_triggered();
    void on_actionRemove_Edge_triggered();

private:
    Edge *findSelectedEdge(bool create);
    Ui::MainWindow *ui;

    std::string filename = "";
    void save(FILE* file);

    QWidget *parent;
    QGraphicsScene *theScene;

    QList<Node* > nodes;
    QList<Edge* > edges;
};

#endif // MAINWINDOW_H
