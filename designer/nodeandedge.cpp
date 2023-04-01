/**
 * @file nodeandedge.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-04-01
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "nodeandedge.h"

#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QGraphicsSceneMouseEvent>
#include <QtGui/QPainter>
#include <QtWidgets/QStyleOption>
#include <QtCore/QtMath>

Edge::Edge(Node *n1, Node *n2, int id): id(id), data{n1, n2} {
    setAcceptedMouseButtons(Qt::NoButton);
    data[0]->edgeList << this;
    data[1]->edgeList << this;
    adjust();
}

void Edge::adjust() {
    if (!data[0] || !data[1])
        return;

    QLineF line(mapFromItem(data[0], 0, 0), mapFromItem(data[1], 0, 0));
    qreal length = line.length();

    prepareGeometryChange();

    if (length > qreal(20.)) {
        QPointF edgeOffset((line.dx() * 10) / length, (line.dy() * 10) / length);
        pointData[0] = line.p1() + edgeOffset;
        pointData[1] = line.p2() - edgeOffset;
    } else {
        pointData[0] = pointData[1] = line.p1();
    }
}

QRectF Edge::boundingRect() const {
    if (!data[0] || !data[1])
        return QRectF();

    return QRectF(pointData[0], QSizeF(pointData[1].x() - pointData[0].x(),
                                    pointData[1].y() - pointData[0].y()))
        .normalized()
        .adjusted(-3, -3, 3, 3);
}

void Edge::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
    if (!data[0] || !data[1])
        return;

    qreal ly = pointData[1].y() - pointData[0].y(), lx = pointData[1].x() - pointData[0].x();
    if (qFuzzyCompare(lx, qreal(0.)) && qFuzzyCompare(ly, qreal(0.)))
        return;

    qreal angle = std::atan2(ly, lx);
    QPointF arrow = QPointF(2 * cos(angle + M_PI / 2), 2 * sin(angle + M_PI / 2));

    painter->setPen(Qt::NoPen);
    painter->setBrush(Qt::darkGray);
    painter->drawPolygon(QPolygonF() << pointData[0] + arrow + QPointF(3, 3) << pointData[0] - arrow + QPointF(3, 3)
                                     << pointData[1] - arrow + QPointF(3, 3) << pointData[1] + arrow + QPointF(3, 3));

    painter->setBrush(Qt::red);
    painter->setPen(QPen(Qt::black, 0));
    painter->drawPolygon(QPolygonF() << pointData[0] + arrow << pointData[0] - arrow
                                     << pointData[1] - arrow << pointData[1] + arrow);
}

Edge::~Edge() {
    for (int i=data[0]->edgeList.size()-1; i>=0; i--)
        if (data[0]->edgeList[i] == this)
            data[0]->edgeList.remove(i);

    for (int i=data[1]->edgeList.size()-1; i>=0; i--)
        if (data[1]->edgeList[i] == this)
            data[1]->edgeList.remove(i);
}

Node::Node(QGraphicsView *graphWidget, int id): id(id), graph(graphWidget) {
    setFlag(ItemIsMovable);
    setFlag(ItemSendsGeometryChanges);
    setCacheMode(DeviceCoordinateCache);
    setZValue(1);
}

void Node::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
    painter->setPen(Qt::NoPen);
    painter->setBrush(Qt::darkGray);
    painter->drawEllipse(-7, -7, 20, 20);

    QRadialGradient gradient(-3, -3, 10);
    if (selected) {
        gradient.setCenter(3, 3);
        gradient.setFocalPoint(3, 3);
        gradient.setColorAt(1, QColor(QColorConstants::Svg::orange).lighter(150));
        gradient.setColorAt(0, QColorConstants::Svg::orange);
    } else {
        gradient.setColorAt(0, QColorConstants::Svg::orange);
        gradient.setColorAt(1, QColor(QColorConstants::Svg::orange).darker(150));
    }
    painter->setBrush(gradient);

    painter->setPen(QPen(Qt::black, 0));
    painter->drawEllipse(-10, -10, 20, 20);
}

QVariant Node::itemChange(GraphicsItemChange change, const QVariant &value) {
    switch (change) {
    case ItemPositionHasChanged:
        for (Edge *edge : std::as_const(edgeList))
        edge->adjust();
        break;
    default:
        break;
    };

    return QGraphicsItem::itemChange(change, value);
}
