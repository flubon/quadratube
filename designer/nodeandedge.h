/**
 * @file nodeandedge.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-04-01
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef NODEANDEDGE_H
#define NODEANDEDGE_H

#include <string>

#include <QtWidgets/QGraphicsItem>
#include <QtWidgets/QGraphicsView>
#include <QtCore/QList>

class Edge;

class Node : public QGraphicsItem {
    friend Edge;
public:
    Node(QGraphicsView *graphWidget, int id);
    inline std::string print() {
        return std::to_string(id) + "\t1\t" + std::to_string(pos().rx()/40) + "\t"
               + std::to_string(pos().ry()/40) + "\t0.0\n";
    }

    enum { Type = UserType + 1 };
    inline int type() const override { return Type; }

    inline QRectF boundingRect() const override {
        qreal adjust = 2;
        return QRectF( -10 - adjust, -10 - adjust, 23 + adjust, 23 + adjust);
    }
    inline QPainterPath shape() const override {
        QPainterPath path;
        path.addEllipse(-10, -10, 20, 20);
        return path;
    }
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;

    bool selected = false;
    QList<Edge *> edgeList;
    int id;

protected:
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override;

    inline void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override {
        selected = !selected;
        update();
        QGraphicsItem::mouseDoubleClickEvent(event);
    }

private:
    QGraphicsView *graph;
};

class Edge : public QGraphicsItem {
    friend Node;
public:
    Edge(Node *n1, Node *n2, int id);
    ~Edge();

    inline Node *operator[](const int i) const { return data[i]; }
    void adjust();

    enum { Type = UserType + 2 };
    inline int type() const override { return Type; }

    inline std::string print() {
        return std::to_string(id) + "\t1\t" + std::to_string(data[0]->id) + "\t" +
               std::to_string(data[1]->id) + "\n";
    }

    inline bool isConnected(Node *node) {
        if (data[0] == node || data[1] == node)
            return true;
        return false;
    }
    int id;

protected:
    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;

private:
    Node *data[2];
    QPointF pointData[2];
};

#endif // NODEANDEDGE_H
